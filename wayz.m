
function [sys,x0,str,ts] = wayz(t,x,u,flag)
switch flag
  case 0
    [sys,x0,str,ts]=mdlInitializeSizes; 
  case 3
    sys=mdlOutputs(t,x,u);
  case 2  % 更新离散状态
    sys = mdlUpdate(t,x,u);
  case {4,9}
    sys=[];
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

function [sys,x0,str,ts]=mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 4;
    sizes.NumOutputs     = 6;
    sizes.NumInputs      = 0;
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 0;   % at least one sample time is needed
    
    sys = simsizes(sizes);
    x0  = [0, 0, 0, 0];
    str = [];
    ts  = [];

function sys = mdlUpdate(t,x,u)
    % 状态更新函数
    prev_t = x(1);    % 前次时间
    
    rost1 = rostime('now');
    t1 = double(rost1.Sec) + double(rost1.Nsec) * 1e-9;
    if(prev_t == 0)
        dt = 0.05;
    else
        dt = t1 - prev_t; %t1 - prev_t; %t - prev_t;
        dt = 0.05;
    end
    T = x(2) + dt;
    
    new_prev_t = t1;

    tadd = x(4) + dt;

    sys = [new_prev_t, T, dt, tadd];

function sys=mdlOutputs(t,x,u)
    
    T = x(2);
    dt = x(3);
    % %定点悬停
    % xk1 = 100;
    % yk1 = -100;
    % zk1 = -30;
    % phi = 0;
    % the = 0;
    % psi = 60; %出角度

    %螺旋线轨迹跟踪
    % % xk1pre = 5*sin(0.2*(T - dt));
    % % yk1pre = 5-5*cos(0.2*(T - dt));
    % xk1 = 25 + 10*sin(0.1*(T));
    % yk1 = 10-10*cos(0.1*(T));
    % psi = 0;%atan2(yk1 - yk1pre, xk1 - xk1pre)/pi*180;
    % zk1 = -1-0.3*(T);
    % phi = 0;
    % the = 0;
    % if (T>150)
    %     cT = 150;
    %     xk1 = 25 + 10*sin(0.1*(cT));
    %     yk1 = 10-10*cos(0.1*(cT));
    %     psi = atan2(0,0);
    %     zk1 = -1-0.3*(cT);
    %     phi = 0;
    %     the = 0;
    % end

    % 直线跟踪
    xk1= 25 + 0.30*T;
    yk1= 0.30*T;
    zk1= -1-0.30*T;   
    phi= 0;
    the= 0;
    psi= 170;

    if(T>190)
        cT = 190; 
        xk1=25 + 0.30*cT;
        yk1=0.30*cT;
        zk1=-1-0.30*120;  
        psi = 170;  
        phi = 0;
        the = 0; 
    elseif (T>120)
        cT = 120; 
        xk1=25 + 0.30*T;
        yk1=0.30*T;
        zk1=-1-0.30*cT;  
        psi = 170;  
        phi = 0;
        the = 0; 
    end
    
    % %for ijfs compare
    % %for test
    % xk1 = 20;
    % yk1 = 0;
    % zk1 = -1.18; % for step response
    % zk1 = -1.35; % for sin curve
    % 
    % if (T>120)
    %     %for test
    %     xk1 = 20;
    %     yk1 = 0;
    % 
    %     psi = 45;
    % 
    %     % zk1 = -1.95; % for step response
    %     zk1 = -2.2 - sin(0.01*2*pi*x(4)-0.2527); % for sin curve
    % end
    
    sys(1) = xk1;
    sys(2) = yk1;
    sys(3) = zk1;
    sys(4) = phi;
    sys(5) = the;
    sys(6) = psi;
