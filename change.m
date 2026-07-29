function [sys,x0,str,ts]=change(t,x,u,flag)
    switch flag
    case 0
        [sys,x0,str,ts]=mdlInitializeSizes;    
    case 3
        sys=mdlOutputs(t,x,u);
    case {2, 4, 9 }
        sys = [];
    otherwise
        error(['Unhandled flag = ',num2str(flag)]);
    end


function [sys,x0,str,ts]=mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 3;
    sizes.NumInputs      = 4;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    
    sys=simsizes(sizes);
    x0=[];
    str=[];
    
    ts=[-1,0];

    
function sys=mdlOutputs(t,x,u)
    if(isnan(u(1)) || (u(1) == 0 && u(2) == 0 && u(3) == 0 && u(4)==0))
        sys=[0 0 0];
    else
        eulerangle = quat2eul([u(1),u(4),u(3),u(2)],'ZYX');
        sys(1) = eulerangle(1) / pi * 180;
        sys(2) = eulerangle(2) / pi * 180;
        sys(3) = eulerangle(3) / pi * 180;
    end
