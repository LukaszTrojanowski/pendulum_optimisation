function [ out ] = Z_cnt(x_in,par)

m = par(1);           % mass of pendulum
M = par(2);           % mass of cart
l = par(3);           % distance of pendulums center of mass from rotation point
b = par(4);           % cart friction coeffitient 
c = par(5);           % pendulum friction coeffitient
J = par(6);           % moment of inertia
g = par(7);           % gravity
k = par(8);

psi_1 = x_in(1);
psi_2 = x_in(2);
psi_3 = x_in(3);
psi_4 = x_in(4);
theta = x_in(5);
theta_dot = x_in(6);
x = x_in(7);
x_dot = x_in(8);
F = x_in(9);

D = l^2*m^2;
A = M + m;
B = D*theta_dot^2;
C = (J*(A) + D*sin(theta)^2 + M*l^2*m);
E = l*m;


psi = [(psi_2*(B*cos(theta)^2 - B*sin(theta)^2 + g*l*cos(theta)*(A) - F*E*sin(theta) + b*E*x_dot*sin(theta)))/C - (psi_4*(E*theta_dot^2*cos(theta) - (E*sin(theta)*(c*theta_dot*(A) + g*l*sin(theta)*(A) + F*E*cos(theta) + B*cos(theta)*sin(theta) - b*E*x_dot*cos(theta)))/C + (E*cos(theta)*(B*cos(theta)^2 - B*sin(theta)^2 + g*l*cos(theta)*(A) - F*E*sin(theta) + b*E*x_dot*sin(theta)))/C - (2*l^3*m^3*cos(theta)^2*sin(theta)*(c*theta_dot*(A) + g*l*sin(theta)*(A) + F*E*cos(theta) + B*cos(theta)*sin(theta) - b*E*x_dot*cos(theta)))/C^2))/(A) - (2*D*psi_2*cos(theta)*sin(theta)*(c*theta_dot*(A) + g*l*sin(theta)*(A) + F*E*cos(theta) + B*cos(theta)*sin(theta) - b*E*x_dot*cos(theta)))/C^2;
       (psi_2*(c*(A) + 2*D*theta_dot*cos(theta)*sin(theta)))/C - (psi_4*(2*E*theta_dot*sin(theta) + (E*cos(theta)*(c*(A) + 2*D*theta_dot*cos(theta)*sin(theta)))/C))/(A) - psi_1;
       0;
       (psi_4*(b + (b*D*cos(theta)^2)/C))/(A) - psi_3 - (b*E*psi_2*cos(theta))/C;
       ];
   
dxdt = [theta_dot;
       (-(A)*c*theta_dot-(A)*g*l*sin(theta)-B*sin(theta)*cos(theta)+m*l*b*x_dot*cos(theta)-m*l*cos(theta)*F)/(J*(A)+m*M*l^2+m^2*l^2*sin(theta)^2);
        x_dot;
       (F - b*x_dot + E*theta_dot^2*sin(theta) + (E*cos(theta)*(c*theta_dot*(A) + g*l*sin(theta)*(A) + F*E*cos(theta) + B*cos(theta)*sin(theta) - b*E*x_dot*cos(theta)))/...
       C)/(A)
       ];
   
z   = (psi_4*((D*cos(theta)^2)/C + 1))/(A) - 2*F*k - (E*psi_2*cos(theta))/C;

out = [psi;dxdt;z];
end

