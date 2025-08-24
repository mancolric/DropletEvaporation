function [u_n, Iteration, flag] = NewtonRaphson(fun, x0, TolF, TolX, MaxIter)

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

        if Iteration==MaxIter
            warning('Reached maximum of iterations')
            flag   = -1;
            break
        end

        Iteration   = Iteration + 1;
        u_n         = u_n + Delta_u_n;

    end
end

