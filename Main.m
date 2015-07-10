%PRZYGOTOWANIE ŒRODOWISKA--------------------------------------------------
clear all;
format long e;
format compact;
%DEFINICJE STA£YCH SYMULACJI-----------------------------------------------
t0 = 0;            %Czas startowy            
tk = 5;            %Czas koñcowy
n = 1000;          %Iloœæ próbek
h = ((tk-t0)/n);   %krok
k = 120;           %Wspó³czynnik funkcji kary
x0 = [0; 0; 0; 0]; %Punkt pocz¹tkowy
u0 = zeros(50, 1); %Sterowanie pocz¹tkowe
u_max = 250;        %Maksymalne dopuszczalne sterowanie
u_min = -250;       %Minimalne dopuszczalne sterowanie
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
%--------------------------------------------------------------------------

%Wykonanie programu optymalizuj¹cego---------------------------------------
[U_in, ~] = BFGS(u0, 10^(-5), 10^(-11), 10^(-4),param,pend_param, u_min,...
    u_max);
[Q,G,psi_in,rk4_out] = rk4b( U_in, param, pend_param);
%--------------------------------------------------------------------------

%Wyniki pracy programu-----------------------------------------------------
figure(1);
Plot(psi_in,rk4_out,U_in,n,3)
figure(2);
Simulator(rk4_out,4,n+1)
%--------------------------------------------------------------------------

