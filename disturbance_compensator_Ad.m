function [sys,x0,str,ts] = disturbance_compensator_Ad(t,x,u,flag,params)
    params.dt    = 0.05;     
    params.u_min = -1500*ones(6,1); 
    params.u_max =  1500*ones(6,1);
    
    params.L_base = diag([2.0, 2.0, 2.0, 0.5, 0.5, 0.5]) / 8;  
    params.L_max  = diag([4.0, 4.0, 4.0, 1.5, 1.5, 1.5]);      
    
    params.gamma  = [0.80; 0.80; 0.80; 0.92; 0.92; 0.92]; 

    switch flag
        case 0
            [sys,x0,str,ts] = mdlInitializeSizes(params);
        case 2
            sys = mdlUpdate(t,x,u,params);
        case 3
            sys = mdlOutputs(t,x,u,params);
        case {1,4,9}
            sys = [];
        otherwise
            error(['Unhandled flag = ', num2str(flag)]);
    end
end

function [sys,x0,str,ts] = mdlInitializeSizes(params)
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 24;    
    sizes.NumOutputs     = 18;    
    sizes.NumInputs      = 27;    
    sizes.DirFeedthrough = 1; 
    sizes.NumSampleTimes = 1;
    
    sys = simsizes(sizes);
    x0  = zeros(24,1); 
    str = [];
    ts  = [-1 0]; 
end

function M_est = get_M_est()
    m = 1863; Xud = -779.79; Yvd = -1222; Zwd = -3659.9;
    Kpd = -534.9; Mqd = -842.69; Nrd = -224.32;
    Ix = 525.39; Iy = 794.20; Iz = 691.23;
    zg = -0.3; 
    
    M_est = zeros(6,6);
    M_est(1,1) = m - Xud; M_est(2,2) = m - Yvd; M_est(3,3) = m - Zwd;
    M_est(4,4) = Ix - Kpd; M_est(5,5) = Iy - Mqd; M_est(6,6) = Iz - Nrd;
    M_est(1,5) = m * zg; M_est(5,1) = m * zg;
    M_est(2,4) = -m * zg; M_est(4,2) = -m * zg;
end

function sys = mdlOutputs(t,x,u,params)
    v_prev       = x(1:6);
    F_nom_prev   = x(7:12);
    d_hat_prev   = x(13:18);
    tau_out_prev = x(19:24);
    
    u_mpc  = u(1:6);      
    v_meas = u(7:12);      
    
    M_est = get_M_est();
    dt = params.dt;
    
    norm_v_lin = norm(v_meas(1:3));
    norm_v_ang = norm(v_meas(4:6));
    
    c_lin = 0.5; 
    c_ang = 0.2; 
    
    % alpha_lin = (norm_v_lin^2) / (norm_v_lin^2 + c_lin^2);
    % alpha_ang = (norm_v_ang^2) / (norm_v_ang^2 + c_ang^2);
    c_vec = [0.3; 0.3; 0.3; 0.15; 0.15; 0.15]; 

    L_k = zeros(6,6);
    % for i = 1:3
    %     L_k(i,i) = params.L_base(i,i) + alpha_lin * (params.L_max(i,i) - params.L_base(i,i));
    % end
    % for i = 4:6
    %     L_k(i,i) = params.L_base(i,i) + alpha_ang * (params.L_max(i,i) - params.L_base(i,i));
    % end
     for i = 1:6
        % 当前通道的速度平方
        v_i_sq = v_meas(i)^2;
        % 当前通道的独立激活系数
        alpha_i = v_i_sq / (v_i_sq + c_vec(i)^2);
        % 计算当前通道的增益
        L_k(i,i) = params.L_base(i,i) + alpha_i * (params.L_max(i,i) - params.L_base(i,i));
     end
     
    d_hat_raw = d_hat_prev + L_k * M_est * (v_meas - v_prev) ...
                - dt * L_k * F_nom_prev - dt * L_k * d_hat_prev;
            
    tau_max =  1000 * ones(6,1); 
    tau_min = -1000 * ones(6,1);
    d_hat_sat = max(min(d_hat_raw, tau_max), tau_min);
    
    tau_out = zeros(6,1);
    for i = 1:6
        tau_out(i) = params.gamma(i) * tau_out_prev(i) + (1 - params.gamma(i)) * d_hat_sat(i);
    end
    
    u_total = u_mpc - tau_out;
    u_total = max(u_total, params.u_min(:));
    u_total = min(u_total, params.u_max(:));
    
    sys = [u_total; tau_out; tau_out ./ diag(M_est)];
end

function sys = mdlUpdate(t,x,u,params)
    v_prev       = x(1:6);
    F_nom_prev   = x(7:12);
    d_hat_prev   = x(13:18);
    tau_out_prev = x(19:24);
    
    u_mpc  = u(1:6);      
    v_meas = u(7:12);      
    euler  = u(13:15);     
    
    M_est = get_M_est();
    dt = params.dt;
    
    norm_v_lin = norm(v_meas(1:3));
    norm_v_ang = norm(v_meas(4:6));
    c_lin = 0.5; 
    c_ang = 0.2; 
    alpha_lin = (norm_v_lin^2) / (norm_v_lin^2 + c_lin^2);
    alpha_ang = (norm_v_ang^2) / (norm_v_ang^2 + c_ang^2);
    
    L_k = zeros(6,6);
    for i = 1:3
        L_k(i,i) = params.L_base(i,i) + alpha_lin * (params.L_max(i,i) - params.L_base(i,i));
    end
    for i = 4:6
        L_k(i,i) = params.L_base(i,i) + alpha_ang * (params.L_max(i,i) - params.L_base(i,i));
    end
    
    d_hat_raw = d_hat_prev + L_k * M_est * (v_meas - v_prev) ...
                - dt * L_k * F_nom_prev - dt * L_k * d_hat_prev;
            
    tau_max =  1000 * ones(6,1); 
    tau_min = -1000 * ones(6,1);
    d_hat_sat = max(min(d_hat_raw, tau_max), tau_min);
    
    tau_out = zeros(6,1);
    for i = 1:6
        tau_out(i) = params.gamma(i) * tau_out_prev(i) + (1 - params.gamma(i)) * d_hat_sat(i);
    end
    
    u_total = u_mpc - tau_out;
    u_total = max(u_total, params.u_min(:));
    u_total = min(u_total, params.u_max(:));
    
    dv_nom = agentdynamic(v_meas, u_total, euler); 
    F_nom_k = M_est * dv_nom;
    
    sys = [v_meas; F_nom_k; d_hat_sat; tau_out]; 
end
