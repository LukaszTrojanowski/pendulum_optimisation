function dxdt = RHS(x, F, par)

m = par(1);           % mass of pendulum
M = par(2);           % mass of cart
l = par(3);           % distance of pendulums center of mass from rotation point
b = par(4);           % cart friction coeffitient 
c = par(5);           % pendulum friction coeffitient
I = par(6);           % moment of inertia
g = par(7);           % gravity
k = par(8);

dxdt = [x(2);
       (-(M+m)*c*x(2)-(M+m)*g*l*sin(x(1))-m^2*l^2*x(2)^2*sin(x(1))*cos(x(1))+m*l*b*x(4)*cos(x(1))-m*l*cos(x(1))*F)/(I*(m+M)+m*M*l^2+m^2*l^2*sin(x(1))^2);
        x(4);
       (F - b*x(4) + l*m*x(2)^2*sin(x(1)) + (l*m*cos(x(1))*(c*x(2)*(M + m) + g*l*sin(x(1))*(M + m) + F*l*m*cos(x(1)) + l^2*m^2*x(2)^2*cos(x(1))*sin(x(1)) - b*l*m*x(4)*cos(x(1))))/...
       (I*(M + m) + l^2*m^2*sin(x(1))^2 + M*l^2*m))/(M + m)
       ];
end
