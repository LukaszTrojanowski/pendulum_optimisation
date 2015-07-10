t0 = 0;            %Czas startowy            
tk = 5;            %Czas koñcowy
n = 1000;          %Iloœæ próbek
h = ((tk-t0)/n);   %krok
k = 120;           %Wspó³czynnik funkcji kary
x0 = [0; 0; 0; 0]; %Punkt pocz¹tkowy
u0 = zeros(50, 1); %Sterowanie pocz¹tkowe
u_max = 80;        %Maksymalne dopuszczalne sterowanie
u_min = -80;       %Minimalne dopuszczalne sterowanie
index = n+1;       
param = [t0 tk n h k index];
%--------------------------------------------------------------------------

%Definicje parametrów wahad³a----------------------------------------------
m = 0.2;           % mass of pendulum
M = 1;             % mass of cart
l = 0.3;           % distance of pendulums center of mass from rotation point
b = 0.1;           % cart friction coeffitient 
c = 0.05;          % pendulum friction coeffitient
J = 0.006;         % moment of inertia
g = 9.8;           % gravity
k = 2;
pend_param = [m M l b c J g k];

[Q, grad, psi_in, rk4_out] = rk4b(zeros(50, 1), param, pend_param);
ep = 1e-6;
dif = [];
for i = 1:length(grad)
   ui_ = zeros(50, 1);
   ui_(i) = ui_(i) + ep;
   [Q_, grad_, ~,~,] = rk4b(ui_, param, pend_param);
   dif(i, 1)= (Q_ - Q)/ep;
end
[dif, grad]