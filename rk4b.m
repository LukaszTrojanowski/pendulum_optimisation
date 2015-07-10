function [Q,G,psi_in,rk4_out] = rk4b( U_in, param, pend_param)

%DEFINICJE STA£YCH SYMULACJI-----------------------------------------------
t0 = param(1);                     
tk = param(2);
n  = param(3);
h  = param(4);
k  = param(5);
index = param(6);
x0 = [0; 0; 0; 0];
%--------------------------------------------------------------------------

%Wektor zmiennych stanu----------------------------------------------------
rk4_out = rk4 ( 'RHS', t0, x0, tk, n, h, U_in, pend_param);
%--------------------------------------------------------------------------

%Warunki pocz¹tkowe rozwi¹zywania pochodnej wskaznika jakoœci--------------
G = zeros(1,50); 
arg_in(:,n+1) = [-2*k*((rk4_out(:,n+1)-[pi;0;0;0]));rk4_out(:,n+1);0];
%--------------------------------------------------------------------------

%Wyznaczanie pochodnych cz¹stkowych wskaznika jakoœci----------------------
for j = 50:-1:1
    arg_in(5:8,index) = rk4_out(:,index);
    arg_in(9,index) = 0;
    for i = 20:-1:1
       k1 =  h * feval('Z_cnt',arg_in(:,index),pend_param);
       arg_tmp = arg_in(:,index) - k1/2;
       k2 =  h * feval('Z_cnt',arg_tmp,pend_param);
       arg_tmp = arg_in(:,index) - k2/2;
       k3 =  h * feval('Z_cnt',arg_tmp,pend_param);
       arg_tmp = arg_in(:,index) - k3;
       k4 =  h * feval('Z_cnt',arg_tmp,pend_param);

       arg_in(:,index-1) = arg_in(:,index) - (k1 + 2*k2 + 2*k3 + k4)/6;
       index = index -1;
    end
    G(j) = arg_in(9,index);
end 
%--------------------------------------------------------------------------

psi_in = arg_in(1:4,:);
W = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];
%Wskaznik jakoœci----------------------------------------------------------
Q = sum(U_in.^2)/1000*5 + k*((rk4_out(:,n+1)-[pi;0;0;0])'*W*(rk4_out(:,n+1)-[pi;0;0;0]));
G = G';
end