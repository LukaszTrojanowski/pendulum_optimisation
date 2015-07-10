function Plot(psi,rk4,F,n,par)

if par == 1
    w  = psi(1,:);
    figure;
    plot(w)
    axis([1 n+1 -inf inf])

    w  = psi(2,:);
    figure;
    plot(w)
    axis([1 n+1 -inf inf])
    hold on;
    
    w  = psi(3,:);
    figure;
    plot(w)
    axis([1 n+1 -inf inf])
    hold on;
    
    w  = psi(4,:);
    figure;
    plot(w) 
    axis([1 n+1 -inf inf])

end

if par == 2
    w  = rk4;
    figure;
    plot(w')
end

if par == 3
    figure;
    u = [];
    for i = 1:length(F)
        u = [u; F(i)*ones(20, 1)];
    end
    plot(u)
end

end

