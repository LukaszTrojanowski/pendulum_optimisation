syms M m l c b F g x x_dot theta theta_dot x_ddot J x_help H psi_1 psi_2 psi_3 psi_4 k psi_1_dot psi_2_dot psi_3_dot psi_4_dot

x_help = F^2;
first(theta, theta_dot, x, x_dot)=theta_dot; 
second(theta, theta_dot, x, x_dot)=(-(M+m)*c*theta_dot-(M+m)*g*l*sin(theta)-m^2*l^2*theta_dot^2*sin(theta)*cos(theta)+m*l*b*x_dot*cos(theta)-m*l*cos(theta)*F)/(J*(m+M)+m*M*l^2+m^2*l^2*sin(theta)^2);
third(theta, theta_dot, x, x_dot)=x_dot;
fourth(theta, theta_dot, x, x_dot)=(F - b*x_dot + l*m*theta_dot^2*sin(theta) + (l*m*cos(theta)*(c*theta_dot*(M + m) + g*l*sin(theta)*(M + m) + F*l*m*cos(theta) + l^2*m^2*theta_dot^2*cos(theta)*sin(theta) - b*l*m*x_dot*cos(theta)))/(J*(M + m) + l^2*m^2*sin(theta)^2 + M*l^2*m))/(M + m);

H(theta, theta_dot, x, x_dot) = psi_1*first(theta, theta_dot, x, x_dot) + psi_2*second(theta, theta_dot, x, x_dot) + psi_3*third(theta, theta_dot, x, x_dot) + psi_4*fourth(theta, theta_dot, x, x_dot) - k*F^2;
%diff(first, theta)
%diff(first, theta_dot)
psi_1_dot = diff(-H, theta);
psi_2_dot = diff(-H, theta_dot);
psi_3_dot = diff(-H, x);
psi_4_dot = diff(-H,x_dot);


z_dot = diff(H,F);

psi_5_dot = 0;
psi_5 = -1;
