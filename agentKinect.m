function dxdt = agentKinect(phi,the,psi,uv)
    u = uv(1);
    v = uv(2);
    w = uv(3);
    p = uv(4);
    q = uv(5);
    r = uv(6);
    Xd=w*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(the)) - v*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(the)) + u*cos(psi)*cos(the);
    Yd=v*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(the)) - w*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(the)) + u*cos(the)*sin(psi);
    Zd=w*cos(phi)*cos(the) - u*sin(the) + v*cos(the)*sin(phi);
    phid=p + r*cos(phi)*tan(the) + q*sin(phi)*tan(the);
    thed=q*cos(phi) - r*sin(phi);
    psid=r*cos(phi)/cos(the) + q*sin(phi)/cos(the);
    dxdt = [Xd; Yd; Zd; phid; thed; psid];
end
