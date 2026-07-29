function [sys,x0,str,ts]=VelBodyChange(t,x,u,flag)
    switch flag
        case 0
            [sys,x0,str,ts]=mdlInitializeSizes;
        case 2
            sys=mdlUpdate(t,x,u);
        case 3
            sys=mdlOutputs(t,x,u);
        case {4, 9 }
            sys = [];
        otherwise
            error(['Unhandled flag = ',num2str(flag)]);
    end

function [sys,x0,str,ts]=mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 6;
    sizes.NumOutputs     = 6;
    sizes.NumInputs      = 10;
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 1;
    
    sys = simsizes(sizes);
    x0 = [0 0 0 0 0 0];
    str = [];
    
    ts = [-1,0];

function sys=mdlUpdate(t,x,u)
    
    sys(1:6) = u(1:6);

    quat = u(7:10)';  %w,x,y,z
    if(norm(quat) == 0)
        quat =[1,0,0,0]; %w,x,y,z
    end
    
    global_vel = sys(1:3);         % 前3项为线速度 [vx_global, vy_global, vz_global]
    global_ang_vel = sys(4:6);     % 后3项为角速度 [wx_global, wy_global, wz_global]
    
    quatforRot = [quat(1), quat(2), quat(3), quat(4)]; %wxyz
    R = quat2rotm(quatforRot);           % 需要 Robotics System Toolbox
    % quatforRot = [quat(2), quat(3), quat(4), quat(1)];  %x,y,z,w  
    % R = quaternion_matrix(quatforRot);
    % R = R(1:3,1:3);

    % 线速度转换：v_body = R^T * v_global
    body_lin_vel = R' * global_vel';
    
    % 角速度转换：ω_body = R^T * ω_global
    body_ang_vel = R' * global_ang_vel';
    
    % 合并结果
    body_vel = [body_lin_vel; body_ang_vel];
    sys = body_vel;

function sys=mdlOutputs(t,x,u)
    sys=x;
