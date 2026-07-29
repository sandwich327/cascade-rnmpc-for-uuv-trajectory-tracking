function dxdt = agentdynamic(vel, uv, euler)
    % vel = [u v w p q r]^T
    u = vel(1);
    v = vel(2);
    w = vel(3);
    p = vel(4);
    q = vel(5);
    r = vel(6);

    % uv = [Tx Ty Tz K M N]^T
    Tx = uv(1);
    Ty = uv(2);
    Tz = uv(3);
    K  = uv(4);
    M  = uv(5);
    N  = uv(6);

    % euler = [phi the psi]^T
    phi = euler(1)*pi/180;
    the = euler(2)*pi/180;
    psi = euler(3)*pi/180; %#ok<NASGU>
    m = 1863;           % kg
    W = m*9.81;       % N
    B = 1.83826*1025*9.81;   % N
    Ix = 525.39;        % kg*m^2
    Iy = 794.20;        % kg*m^2
    Iz = 691.23;        % kg*m^2

    rg = [0, 0, 0];    % m
    rb = [0, 0, -0.3];  % m
    BG = rg - rb;      % m
    xg = BG(1);        % m
    yg = BG(2);        % m
    zg = BG(3);        % m

    Xud = -779.79;   % kg
    Yvd = -1222;     % kg
    Zwd = -3659.9;   % kg
    Kpd = -534.9;    % kg*m^2/rad
    Mqd = -842.69;   % kg*m^2/rad
    Nrd = -224.32;   % kg*m^2/rad

    Xu = -74.82;   % N*s/m
    Yv = -69.48;   % N*s/m
    Zw = -782.4;   % N*s/m
    Kp = -268.8;   % N*s/rad
    Mq = -309.77;   % N*s/rad
    Nr = -105;      % N*s/rad

    Xuu = -748.22;  % N*s^2/m^2
    Yvv = -992.53;  % N*s^2/m^2
    Zww = -1821.01; % N*s^2/m^2
    Kpp = -672;     % N*s^2/rad^2
    Mqq = -774.44;  % N*s^2/rad^2
    Nrr = -523.27;  % N*s^2/rad^2

    ud=(m*zg*(M - r*(Ix*p - Kpd*p) + p*(Iz*r - Nrd*r) + w*(Xud*u - m*u) - u*(Zwd*w - m*w) + q*(Mq + Mqq*abs(q)) - W*zg*sin(the)))/(Iy*Xud - Mqd*Xud + m^2*zg^2 - Iy*m + Mqd*m) - ((Iy - Mqd)*(Tx + u*(Xu + Xuu*abs(u)) + sin(the)*(W - B) - r*(Yvd*v - m*v) + q*(Zwd*w - m*w)))/(Iy*Xud - Mqd*Xud + m^2*zg^2 - Iy*m + Mqd*m);
    vd=-((Ix - Kpd)*(Ty + v*(Yv + Yvv*abs(v)) + r*(Xud*u - m*u) - p*(Zwd*w - m*w) - cos(the)*sin(phi)*(W - B)))/(Ix*Yvd - Kpd*Yvd + m^2*zg^2 - Ix*m + Kpd*m) - (m*zg*(K + r*(Iy*q - Mqd*q) - q*(Iz*r - Nrd*r) - w*(Yvd*v - m*v) + v*(Zwd*w - m*w) + p*(Kp + Kpp*abs(p)) - W*zg*cos(the)*sin(phi)))/(Ix*Yvd - Kpd*Yvd + m^2*zg^2 - Ix*m + Kpd*m);
    wd=-(Tz + w*(Zw + Zww*abs(w)) - q*(Xud*u - m*u) + p*(Yvd*v - m*v) - cos(phi)*cos(the)*(W - B))/(Zwd - m);
    pd=((Yvd - m)*(K + r*(Iy*q - Mqd*q) - q*(Iz*r - Nrd*r) - w*(Yvd*v - m*v) + v*(Zwd*w - m*w) + p*(Kp + Kpp*abs(p)) - W*zg*cos(the)*sin(phi)))/(Ix*Yvd - Kpd*Yvd + m^2*zg^2 - Ix*m + Kpd*m) - (m*zg*(Ty + v*(Yv + Yvv*abs(v)) + r*(Xud*u - m*u) - p*(Zwd*w - m*w) - cos(the)*sin(phi)*(W - B)))/(Ix*Yvd - Kpd*Yvd + m^2*zg^2 - Ix*m + Kpd*m);
    qd=((Xud - m)*(M - r*(Ix*p - Kpd*p) + p*(Iz*r - Nrd*r) + w*(Xud*u - m*u) - u*(Zwd*w - m*w) + q*(Mq + Mqq*abs(q)) - W*zg*sin(the)))/(Iy*Xud - Mqd*Xud + m^2*zg^2 - Iy*m + Mqd*m) + (m*zg*(Tx + u*(Xu + Xuu*abs(u)) + sin(the)*(W - B) - r*(Yvd*v - m*v) + q*(Zwd*w - m*w)))/(Iy*Xud - Mqd*Xud + m^2*zg^2 - Iy*m + Mqd*m);
    rd=(N + r*(Nr + Nrr*abs(r)) + q*(Ix*p - Kpd*p) - p*(Iy*q - Mqd*q) - v*(Xud*u - m*u) + u*(Yvd*v - m*v))/(Iz - Nrd);
    dxdt = [ud; vd; wd; pd; qd; rd];
end
