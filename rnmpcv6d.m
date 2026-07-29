function [sys,x0,str,ts] = rnmpcv6d(t,x,u,flag,fuzzyOn)
    persistent t_calc_history
    switch flag
        case 0
            [sys,x0,str,ts] = mdlInitializeSizes;
            t_calc_history = [];
        case 3  % 输出（计算控制量）
            [sys, t_calc_history] = mdlOutputs(t,x,u,fuzzyOn, t_calc_history);
        case {2,4}
            sys = [];
        case 9
            assignin('base', 't_calc_history', t_calc_history);
        otherwise
            error(['Unhandled flag = ', num2str(flag)]);
    end
end

%% 初始化
function [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 0;      
    sizes.NumOutputs     = 12;     % 6 维体坐标速度 + 6 error
    sizes.NumInputs      = 18;     % [期望位姿6; 当前位姿6; 当前速度6]
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;

    sys = simsizes(sizes);
    x0  = [];
    str = [];

    ts  = [-1 0];
    simStateCompliance = 'UnknownSimState';
end

%% 输出：求解 NMPC
function [sys, t_calc_history] = mdlOutputs(t, x, u, fuzzyOn, t_calc_history)
    persistent DV_last;

    % u1-6: 期望位姿
    % u7-12: 当前位姿
    % u13-18: 当前体坐标速度

    cN = 2;      % 控制段数
    pN = 8;      % 预测步数

    % 角度转弧度
    u(4:6)   = u(4:6)   * pi / 180;
    u(10:12) = u(10:12) * pi / 180;

    dt = 0.05;  % 20 Hz

    % % 特殊逻辑：当水平距离较远时，强制调整航向指向目标
    % poserr = u(1:2) - u(7:8);
    % if norm(poserr) > 5
    %     heading = atan2(poserr(2), poserr(1));
    %     u(6) = heading;
    % end

    % ---- 1. 计算各通道误差 ----
    pos_err  = u(1:3)  - u(7:9);
    ang_err  = u(4:6)  - u(10:12);
    ang_err  = mod(ang_err + pi, 2*pi) - pi; % 角度归一化到 [-pi, pi]
    
    % 组合成 6x1 的误差向量
    error_vec = [pos_err; ang_err];

    % ---- 2. 基准 Q、R 对角线元素 ----
    Q_base_diag = [4 4 5 14 14 18];
    R_base_diag = [1 1 1 2 2 2] * 3.5/2;
    % R_base_diag(3:6) = R_base_diag(3:6) * 0.9;
    % ---- 3. 独立模糊调整 (核心修改) ----
    % 初始化缩放系数向量
    alphaQ_vec = zeros(1, 6);
    alphaR_vec = zeros(1, 6);

    % 遍历 6 个自由度，分别计算权重
    for i = 1:6        
        if(fuzzyOn == 1)
            % 传入该通道的绝对误差
            abs_err_i = abs(error_vec(i));
            if i>=3
                abs_err_i = 3*abs_err_i;
            end
            
            % 调用模糊函数
            [aQ, aR] = fuzzyQR_gain(abs_err_i);
            
            alphaQ_vec(i) = aQ;
            alphaR_vec(i) = aR;
        else
            alphaQ_vec(i) = 1;
            alphaR_vec(i) = 1;
        end
    end

    % 构造最终的对角矩阵
    % Q_ii = alphaQ_i * Q_base_ii
    Q = diag(alphaQ_vec .* Q_base_diag);
    R = diag(alphaR_vec .* R_base_diag);
    
    % 终端代价 F 通常与 Q 保持比例关系
    F = 6 * Q;

    % ---- 后续 NMPC 求解逻辑保持不变 ----
    current_v = u(13:18);   % 6×1

    % 速度边界
    vmin = -1 * ones(6,1);
    vmax =  1 * ones(6,1);

    n_dec = 6 * cN; 

    A = []; b = []; Aeq = []; beq = [];
    dv_max_vec = [0.5 0.5 0.5 0.5 0.5 0.5]';%[0.25 0.25 0.2  0.25 0.2 0.2]';
    lb_v = max(-dv_max_vec, vmin - current_v);
    ub_v = min( dv_max_vec, vmax - current_v);

    lb = repmat(lb_v, cN, 1);
    ub = repmat(ub_v, cN, 1);

    x_current = u(7:12);
    desired_x = u(1:6);

    t_calc_start = tic;

    cost_fun = @(DV) costMPC_2segment(dt, x_current, DV, cN, pN, ...
                                      Q, F, R, desired_x, current_v);

    nonlcon = [];

    options = optimoptions(@fmincon, ...
        'Algorithm', 'sqp', ...
        'Display', 'off', ...
        'MaxIterations', 50, ...
        'SpecifyObjectiveGradient', false, ...
        'SpecifyConstraintGradient', false, ...
        'MaxFunctionEvaluations', 100, ...
        'OptimalityTolerance', 1e-4, ...
        'StepTolerance',       1e-4, ...
        'ConstraintTolerance', 1e-4);

    % warm start
    if isempty(DV_last) || numel(DV_last) ~= n_dec
        DV0 = zeros(n_dec,1);
    else
        u_last = reshape(DV_last, 6, cN);
        u_new  = [u_last(:,2), u_last(:,2)];
        DV0    = u_new(:);
    end
    
    DV_opt = fmincon(cost_fun, DV0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    
    t_calc = toc(t_calc_start);
    if t_calc > dt/2
        fprintf('Warning: MPC solve time = %.3f s\n', t_calc);
    end
    t_calc_history(end+1,1) = t_calc;

    DV_last = DV_opt;

    dv1   = DV_opt(1:6);
    v_out = current_v + dv1;
    v_out = min(max(v_out, vmin), vmax);

    sys(1:6) = v_out;
    sys(7:9) = pos_err;
    sys(10:12) = ang_err;
end

%% 模糊增益函数（运动学 NMPC 用，单通道误差 -> alphaQ, alphaR）
function [alphaQ, alphaR] = fuzzyQR_gain(err_val)
    % 输入 err_val: 单个通道的绝对误差 (标量)，外面已对角度做了缩放
    e = max(err_val, 0);

    %========================
    % 1. 隶属度函数重新设计
    %========================
    % 这里取一个比较保守、平滑的划分：
    % small :   0   ~ 0.4 左右
    % mid   :   0.2 ~ 1.2
    % large :   0.8 ~ 2.5+
    %
    % 你现在的调用里，位置就是米级，角度外面乘了 3 以后大概 0~几的量级，
    % 这个范围可基本覆盖绝大多数情况，如果你实测误差明显更大/更小，
    % 可以在下面微调阈值。

    % --- small，左肩型隶属：0~0.4 ---
    if e <= 0.2
        mu_small = 1;
    elseif e >= 0.6
        mu_small = 0;
    else
        mu_small = (0.6 - e) / (0.6 - 0.2);  % 0.2~0.6 线性下降
    end

    % --- mid，三角形：0.2~1.2 ---
    if e <= 0.2 || e >= 1.2
        mu_mid = 0;
    elseif e <= 0.7
        mu_mid = (e - 0.2) / (0.7 - 0.2);   % 0.2~0.7 上升
    else
        mu_mid = (1.2 - e) / (1.2 - 0.7);   % 0.7~1.2 下降
    end

    % --- large，右肩型：0.8~2.5 ---
    if e <= 0.8
        mu_large = 0;
    elseif e >= 2.5
        mu_large = 1;
    else
        mu_large = (e - 0.8) / (2.5 - 0.8); % 0.8~2.5 线性上升
    end

    % 归一化，避免数值问题
    mu_sum = mu_small + mu_mid + mu_large;
    if mu_sum < 1e-6
        % 极端情况下直接当 small 处理
        mu_small = 1; mu_mid = 0; mu_large = 0;
        mu_sum   = 1;
    end
    mu_small = mu_small / mu_sum;
    mu_mid   = mu_mid   / mu_sum;
    mu_large = mu_large / mu_sum;

    %========================
    % 2. 规则库：Q、R 的趋势（针对扰动/稳定性）
    %========================
    % 设计目标：
    % - 小误差：重视稳态精度，但不过度减小 R，避免高频小抖
    % - 中误差：适度加大 Q，同时稍微降低 R，加快收敛
    % - 大误差：Q 不再继续猛增，只是略大一点；R 不要极小，否则会太冲
    %
    % 注意：这里是“运动学层”Q/R，最终作用是影响预测轨迹形状，
    % 动力学那层你已经单独有一套更保守的模糊逻辑，两层叠加后整体不会太激进。

    % small 区（误差已经较小，主要压扰动、减小稳态偏差）
    alphaQ_small  = 1.2;   % 比 1 稍大一点，稳态段更紧
    alphaR_small  = 1.1;   % R 稍大，抑制小范围乱晃

    % mid 区（典型跟踪过程，允许多动一点）
    alphaQ_mid    = 1.35;  % 合理增大 Q，加快误差收敛
    alphaR_mid    = 0.95;  % R 稍低，允许更积极一点的修正

    % large 区（初始偏差大或外界强扰，既要追又不能太猛）
    alphaQ_large  = 1.5;   % 不再像原先那样继续拉得很高，只比 mid 稍高
    alphaR_large  = 0.85;  % 比 mid 再小一点，但远没到“非常小”的程度

    % 如果你发现“误差大时轨迹太激进”，可以适当调大 alphaR_large，比如 0.85 -> 1.0
    % 如果你觉得“收敛偏慢”，可以略微增大 alphaQ_mid / alphaQ_large

    %========================
    % 3. 解模糊
    %========================
    alphaQ = mu_small * alphaQ_small + ...
             mu_mid   * alphaQ_mid   + ...
             mu_large * alphaQ_large;

    alphaR = mu_small * alphaR_small + ...
             mu_mid   * alphaR_mid   + ...
             mu_large * alphaR_large;
end

%% 代价函数：两段常值控制（DV = [dv1; dv2]）
function cost = costMPC_2segment(dt, x0, DV, cN, pN, Q, F, R, desiredx, v0)
    % x0: 当前状态（6×1）
    % DV: (6*cN)×1 = 12×1
    % Q, F, R: 现在是根据各通道误差动态调整后的对角矩阵

    if cN ~= 2
        error('costMPC_2segment assumes cN = 2.');
    end
    
    k1 = 4;
    k1 = min(k1, pN); 

    dv_mat = reshape(DV, 6, cN); 
    dv1 = dv_mat(:,1);
    dv2 = dv_mat(:,2);

    e = zeros(12*pN, 1); 

    x = x0;
    f = @agentKinect; % 确保你的工作区或路径中有这个函数
    RKN = 1;

    for i = 1:pN
        if i <= k1
            dv_i = dv1;
        else
            dv_i = dv2;
        end

        u = v0 + dv_i;

        % 状态预测
        for j = 1:RKN
            phi = x(4); the = x(5); psi = x(6);
            f1  = f(phi, the, psi, u);
            x   = x + f1 * (dt / RKN);
        end

        % 状态误差
        nuralact      = zeros(6,1);
        nuralact(1:3) = desiredx(1:3) - x(1:3);
        nuralact(4:6) = desiredx(4:6) - x(4:6);
        nuralact(4:6) = mod(nuralact(4:6) + pi, 2*pi) - pi;

        % 计算代价项
        % 注意：这里 Q 和 R 已经是调整好的对角矩阵，直接乘即可
        if i == pN
            % 终端代价 a = F * error
            a = F * nuralact;
        else
            % 过程代价 a = Q * error
            a = Q * nuralact;
        end

        % 控制增量代价 b = R * dv
        b = R * dv_i;

        e((i-1)*12+1 : i*12) = [a; b];
    end

    cost = e' * e;
end
