function  Simulator(x0,par,n)

h_cart = plot(NaN, NaN, 'Marker', 'square', 'color', 'red', 'LineWidth', 6);
hold on
h_pend = plot(NaN, NaN, 'bo', 'LineWidth', 3);
axis([-5 5 -5 5]);
axis manual;
xlim([-25 25]);
ylim([-25 25]);

for i = 1:n
    l = par;
    set(h_cart, 'XData', x0(3,i), 'YData', 0, 'LineWidth', 5);
    set(h_pend, 'XData', sin(x0(1,i))*l+x0(3,i), 'YData', -cos(x0(1,i))*l, 'LineWidth', 2);
    pause(0.003);
end;



end

