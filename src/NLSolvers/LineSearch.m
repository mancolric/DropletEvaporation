function [u_np1, f_np1, fnorm_np1, flagLS] = LineSearch(fun, u_n, Delta_u_n, fnorm_n, LS_max_iteration) 

    %flagLS = -1:   Iteration > Max Iteration

    lambda          = 1;
    LS_iteration    = 0;
    flagLS          = 0;

    while true

        LS_iteration    = LS_iteration + 1;

        if LS_iteration>LS_max_iteration
            warning('Reached maximum of line search iterations')
            flagLS = -1;
            break
        end

        u_np1           = u_n + lambda*Delta_u_n;
        [f_np1, ~]      = fun(u_np1, false);
        fnorm_np1       = norm(f_np1);

        if fnorm_np1<fnorm_n
            break
        else
            lambda      = lambda/2;
        end
        
    end
    
end

