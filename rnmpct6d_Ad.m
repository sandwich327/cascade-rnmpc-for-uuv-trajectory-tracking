function [sys,x0,str,ts] = rnmpct6d_Ad(t,x,u,flag,fuzzyOn)
    persistent tt_calc_history
    switch flag
    case 0
        [sys,x0,str,ts] = mdlInitializeSizes;
        tt_calc_history = [];
    case 3 % 输出（计算控制量）
        [sys, tt_calc_history] = mdlOutputs(t,x,u,fuzzyOn, tt_calc_history);
    case {2,4}
        sys = [];
    case 9
        assignin('base', 'tt_calc_history', tt_calc_history);
    otherwise
        error(['Unhandled flag = ', num2str(flag)]);
    end
end

%% 初始化
function [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 18;  % 6 维体力/力矩 [Tx Ty Tz K M N]
    sizes.NumInputs      = 24; % [期望速度6; 当前速度6; 当前位姿6]
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;

    sys = simsizes(sizes);
    x0  = [];
    str = [];
    ts  = [-1 0];
    simStateCompliance = 'UnknownSimState';
end

%% 输出：求解 NMPC
function [sys, tt_calc_history] = mdlOutputs(t, x, u, fuzzyOn, tt_calc_history)
    persistent DV_last u_out_last x_nom_last K_fb Delta_u e_max_last

    if isempty(K_fb)
        K_fb = [ ...
            11.6162   -0.0000    0.0000   -0.0000   -1.3712   -0.0000;
             0.0000   16.1009    0.0000    1.8758    0.0000    0.0000;
            -0.0000    0.0000   10.5956    0.0000   -0.0000   -0.0000;
            -0.0000   -0.2963    0.0000    2.3114    0.0000   -0.0000;
             0.1902   -0.0000   -0.0000   -0.0000    2.0091    0.0000;
             0.0000   -0.0000   -0.0000   -0.0000   -0.0000   13.9476];

        e_max = [0.25; 0.25; 0.2; 0.2; 0.2; 0.2]*4;

        % 收缩量 Δu = |K_fb| e_max
        Delta_u = abs(K_fb) * e_max;
    end

    % ==== 输入拆分 ====
    vel_des   = u(1:6);        % 期望速度
    vel_cur   = u(7:12);       % 当前速度
    pose_cur  = u(13:18);      % 当前位姿
    euler_cur = pose_cur(4:6); % [phi the psi]    
    a_dist   = u(19:24);      

    e_max_min = [0.25; 0.25; 0.2; 0.2; 0.2; 0.2] * 0.25;
    e_max_max = [0.25; 0.25; 0.2; 0.2; 0.2; 0.2] * 1;
    c_lambda = [0.005; 0.005; 0.003; 0.003; 0.003; 0.003] * 1; 

    w_bar = abs(a_dist);
    e_max_raw = e_max_min + c_lambda .* w_bar;
    e_max_raw = min(max(e_max_raw, e_max_min), e_max_max); 

    if isempty(e_max_last)
        e_max_last = e_max_min;
    end
    alpha_emax = 0.2; % 滤波系数，可调 (0.1~0.5)
    e_max_k = (1 - alpha_emax) * e_max_last + alpha_emax * e_max_raw;
    e_max_last = e_max_k;

    Delta_u = abs(K_fb) * e_max_k;

    if isempty(x_nom_last)
        x_nom_last = vel_cur;         % 初始名义速度 = 当前真实
        u_out_last = zeros(6,1);      % 记录上一次输出控制
        DV_last    = zeros(12,1);     % 6*cN, cN = 2
    end

    % ===== MPC 结构参数 =====
    cN = 2;        % 控制段数
    pN = 6;        % 预测步数
    dt = 0.05;     % 固定采样时间

    % ---- 1. 计算各通道速度误差 ----
    vel_err        = vel_des - vel_cur;
    err_vec_scaled = vel_err * 2; 
        
    err_norm = norm(vel_err);
    alpha_max = 0.95;  
    alpha_min = 0.5;   
    err_ref   = 0.3;   
    s = min(err_norm / err_ref, 1);   
    alpha_smooth = alpha_max - (alpha_max - alpha_min) * s; %#ok<NASGU>

    Q_base_diag = [20 20 40 20 20 40];
    R_base_diag = [1 1 1 2 2 1] * 0.0045/2;

    alphaQ_vec = zeros(1, 6);
    alphaR_vec = zeros(1, 6);

    for i = 1:6
        if fuzzyOn == 1
            abs_err_i = 5 * abs(err_vec_scaled(i));
            [aQ, aR] = fuzzyQR_gain_channel(abs_err_i, i);
            alphaQ_vec(i) = aQ;
            alphaR_vec(i) = aR;
        else
            alphaQ_vec(i) = 1;
            alphaR_vec(i) = 1;
        end
    end

    Q = diag(alphaQ_vec .* Q_base_diag);
    R = diag(alphaR_vec .* R_base_diag);

    F = 4 * Q;     

    n_dec = 6 * cN; 

    u_phys_min = -1500 * ones(6,1);
    u_phys_max =  1500 * ones(6,1);

    umin = u_phys_min + Delta_u;
    umax = u_phys_max - Delta_u;

    A   = []; 
    b   = []; 
    Aeq = []; 
    beq = [];
    lb  = repmat(umin, cN, 1);   
    ub  = repmat(umax, cN, 1);

    x_current = vel_cur;   
    desired_x = vel_des;   

    cost_fun = @(DV) costMPC_2segment( ...
        dt, x_current, DV, cN, pN, Q, F, R, desired_x, euler_cur, K_fb);

    nonlcon = []; 

    options = optimoptions(@fmincon, ...
        'Algorithm', 'sqp', ...
        'Display', 'off', ...
        'MaxIterations', 50, ...
        'MaxFunctionEvaluations', 200, ...
        'OptimalityTolerance', 1e-4, ...
        'StepTolerance',       1e-5, ...
        'ConstraintTolerance', 1e-4);

    % ==== 初始值（warm start）====
    if isempty(DV_last) || numel(DV_last) ~= n_dec
        DV0 = zeros(n_dec,1);
    else
        u_last = reshape(DV_last, 6, cN); 
        u_new  = [u_last(:,2), u_last(:,2)];
        DV0    = u_new(:);
    end

    % ==== 求解 MPC ====
    t_calc_start = tic;
    DV_opt = fmincon(cost_fun, DV0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    t_calc = toc(t_calc_start);

    if t_calc > dt/2
        fprintf('Warning: MPC t = %.3f s\n', t_calc);
    end
    tt_calc_history(end+1,1) = t_calc;

    DV_last = DV_opt;

    % ==== 原始 MPC 输出（第一段控制量，名义控制） ====
    u1    = DV_opt(1:6);                
    u_mpc = min(max(u1, umin), umax);   % 名义控制 u_nom

    % --------- 名义状态滚动 & Tube 反馈 ---------
    f   = @agentdynamic;   
    RKN = 1;               

    x_nom = x_nom_last; %vel_cur; %
    for j = 1:RKN
        dx_nom = f(x_nom, u_mpc, euler_cur);
        x_nom  = x_nom + dx_nom * (dt / RKN);
    end
    x_nom_now = x_nom;

    % tube 反馈：用名义-真实速度误差修正控制
    e_tube   = x_nom_last - vel_cur;  
    u_tube   = u_mpc + K_fb * e_tube;

    % 记录新的名义状态作为下一步初值
    x_nom_last = x_nom_now;

    % ==== 跨时间步的一阶输出平滑 ====
    if isempty(u_out_last)
        u_out_last = u_tube;
    end

    % 当前实现：直接使用 u_tube，可按需要加入 alpha_smooth 平滑
    % u_out = u_out_last + alpha_smooth * (u_tube - u_out_last);
    u_out = u_tube;

    % 保持在收缩约束内
    u_out = min(max(u_out, umin), umax);

    u_out_last = u_out;

    sys(1:6) = u_out;
    sys(7:12) = e_tube;
    sys(13:18) = e_max_k;
end

function cost = costMPC_2segment(dt, x0, DV, cN, pN, Q, F, R, desiredx, euler, K_fb)

    if cN ~= 2
        error('costMPC_2segment assumes cN = 2.');
    end

    % 前 k1 步用 u1，后面用 u2
    k1 = 3; 
    k1 = min(k1, pN);

    u_mat = reshape(DV, 6, cN); 
    u1 = u_mat(:,1);
    u2 = u_mat(:,2);

    % 每步存 [Q*err; R*u] 共 12 维
    e = zeros(12*pN, 1);

    x   = x0;    
    f   = @agentdynamic; 
    RKN = 1;             

    for i = 1:pN
        if i <= k1
            ui = u1;
        else
            ui = u2;
        end

        % === 状态预测 ===
        for j = 1:RKN
            dx  = f(x, ui, euler);
            x   = x + dx * (dt / RKN);
        end

        % === 速度误差 ===
        nuralact = desiredx - x;  

        % 误差加权（末步仍然使用 F）
        if i == pN
            a = F * nuralact;
        else
            a = Q * nuralact;
        end

        % 控制加权
        ui_cost = ui;
        ui_cost(3) = ui_cost(3) + 260; % 深度补偿浮力
        b = R * ui_cost;

        e((i-1)*12+1 : i*12) = [a; b];
    end

    lambda_du    = 1;    
    du12         = u2 - u1;
    smooth_cost  = lambda_du * (du12.' * du12);

    xN = x;                      % 预测末端状态
    eN = desiredx - xN;          % 末端误差（按你原来符号）

    lambda_xT = 1.5;             % 终端状态权重系数，可 0.5~5 调整
    P_term    = lambda_xT * F;   % P ≈ F

    u_fbN = -K_fb * eN;

    lambda_uT = 0.5;             % 可在 0.01~0.5 内尝试
    R_fb      = lambda_uT * R;   % 用同一结构的 R 作缩放

    J_term = eN.' * P_term * eN + u_fbN.' * R_fb * u_fbN;

    cost = e' * e + smooth_cost + J_term;
end
