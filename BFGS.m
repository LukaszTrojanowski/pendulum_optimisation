function [ x, Q ]=BFGS( x0, eps0, eps1, eps2, param, pend_param, u_min,...
    u_max)
%INICJALIZACJA-------------------------------------------------------------
go_flag = 1;
old_x = 0;
old_g = 0;
Q = inf;
x = x0;
[Q,grad] = feval('rk4b', x0,param , pend_param);
disp([Q, norm(grad)]);
%PÊTLA G£ÓWNA--------------------------------------------------------------
while(norm(grad) > eps0)
    
    if (go_flag)
        W = eye(length(x0));
    else 
        r = grad - old_g;
        s = x - old_x;
        if (norm(r) > eps1)
            W = W + (r*r')/(s'*r) - (W*s*s'*W)/(s'*W*s);
        else
            W = eye(length(x0));
        end
    end
    d = -W\grad;
    
    if (d'*grad >= 0)
        W = eye(length(x0));
        [~, grad] = constr_feval('rk4b', x0,param , pend_param, u_min, u_max);
        d = -grad;
    end
    
    old_x = x;
    old_g = grad;
    old_Q = Q;
    [x, Q] = step(x, d, old_Q,param , pend_param, u_min, u_max);

    if (norm(x-old_x) > eps2)
        go_flag = 0;
        [Q10,grad] = constr_feval('rk4b', x,param , pend_param, u_min, u_max);
        disp([Q10, norm(grad), grad']);
    elseif (go_flag == 0)
        go_flag = 1;
    else
        break;
    end
        
end


end

%FUNKCJA WYKONUJ¥CA KROK W KIERUNKU d--------------------------------------
function [x, q] = step(m_x, m_g,  m_q,param , pend_param, u_min, u_max)
    alpha = 1;
    x = m_x;
    q = m_q;
    while(alpha> 1e-24)
        x_temp = m_x + alpha*m_g;
        %SPRAWDZENIE OGRANICZEÑ--------------------------------------------
        constraints = (x_temp > u_max) | (x_temp < u_min);
        if any(constraints)
           x_temp(x_temp > u_max) = u_max;
           x_temp(x_temp < u_min) = u_min;
        end
        %------------------------------------------------------------------
        [temp_q, temp_g] = feval('rk4b', x_temp,param , pend_param);
        %temp_g(constraints) = 0;
        if (temp_q < m_q)
            x = x_temp;
            q = temp_q;
            %g = temp_g;
            break;
        else
            alpha = alpha*0.5;
        end
    end
end
%--------------------------------------------------------------------------

%FUNKCJA OBLICZAJACA GRADIENT Z UWZGLÊDNIENIEM OGRANICZEÑ------------------
function [q, g] = constr_feval(rk, m_x,m_param , m_pend_param, u_min, u_max)
    [temp_q, temp_g] = feval(rk, m_x,m_param , m_pend_param);
    temp_g((m_x >= u_max) | (m_x <= u_min)) = 0;
    q = temp_q;
    g = temp_g;
end
%--------------------------------------------------------------------------