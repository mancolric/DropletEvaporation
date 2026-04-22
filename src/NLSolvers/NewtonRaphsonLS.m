function [u_n, Iteration, flag] = NewtonRaphsonLS(fun, x0, TolF, TolX, MaxIter, LS_max_iteration)

%First iteration:
Iteration   = 0;
u_n         = x0;

%Iterate:
while true

    %This step could be executed only if mod(Iteration, some number)=0:
    [f_n,J_n]   = fun(u_n, true);
    fnorm_n     = norm(f_n);
    
    if fnorm_n<TolF
        flag    = 1;
        return
    end

    Delta_u_n       = - (J_n\f_n);
    if norm(Delta_u_n)<TolX
        flag    = 2;
        break
    end

    if any(isnan(Delta_u_n))
        warning('Found NaNs in preconditioned residual')
        flag    = -3;
        break
    end
    if Iteration==MaxIter
        warning('Reached maximum of iterations')
        flag   = -1;
        break
    end
    
    [u_np1, f_np1, fnorm_np1, flagLS] = LineSearch(fun, u_n, Delta_u_n, fnorm_n, LS_max_iteration); 
    if flagLS<0
        warning('Line search failed')
        flag    = -2;
        break
    end

    Iteration   = Iteration + 1;
    u_n         = u_np1;
    f_n         = f_np1;
    fnorm_n     = fnorm_np1;

end
end

