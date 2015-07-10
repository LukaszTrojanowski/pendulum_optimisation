function [wi, ti] = rk4 ( RHS, t0, x0, tf, n, h, U_in, par)
neqn = length ( x0 );
ti = linspace ( t0, tf, n+1 );
wi = [ zeros( neqn, n+1 ) ];
wi(1:neqn, 1) = x0';

u = [];
for i = 1:length(U_in)
    u = [u; U_in(i)*ones(20, 1)];
end
u = [u;U_in(50)];
for i = 1:n
    k1 = h * feval ( 'RHS', x0, u(i),par);
    k2 = h * feval ( 'RHS', x0 + (k1/2), u(i),par);
    k3 = h * feval ( 'RHS', x0 + k2/2, u(i),par);
    k4 = h * feval ( 'RHS', x0 + k3, u(i),par);
    x0 = x0 + ( k1 + 2*k2 + 2*k3 + k4 ) / 6;
    t0 = t0 + h;
    wi(1:neqn,i+1) = x0';  
end;




