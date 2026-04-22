%[x_n, nIters, flag] = FixedPointIter(fun, x0, TolG, MaxIter)
%
%Normally, an equation of the form
%   f(x)    = 0
%can be preconditioned with an approximate Jacobian J in such a way that 
%it is equivalent to
%   g(x)    = J(x)\f(x) = 0.
%This function solves the equation
%   g(x)    = 0
%via the fixed point iteration
%   x       = x - g(x)
%fun receives the argument x and returns g(x) 
function [x_n, Iteration, flag] = FixedPointIter(fun, x0, TolG, MaxIter)

    %First iteration:
    Iteration   = 0;
    x_n         = x0;
    
    %Iterate:
    while true
        %Preconditioned residual:
        g_n         = fun(x_n);

        %Exit loop:
        if norm(g_n)<=TolG
            flag    = 2;
            break
        elseif Iteration==MaxIter
            warning('Reached maximum of iterations')
            flag   = -1;
            break
        elseif any(isnan(g_n))
            warning('Found NaNs in preconditioned residual')
            flag    = -2;
            break
        end

        %Update solution:
        Iteration   = Iteration + 1;
        x_n         = x_n - g_n;
        


    end
    

    %Display info:
    if true
        if flag>0
            disp(['FixedPoint converged', ...
                ', nIters=', sprintf('%d', Iteration), ...
                ', |x|=', sprintf('%.3E', norm(x_n)), ...
                ', |g|=', sprintf('%.3E', norm(g_n))])
        else
            disp(['FixedPoint failed', ...
                ', nIters=', sprintf('%d', Iteration), ...
                ', |x|=', sprintf('%.3E', norm(x_n)), ...
                ', |g|=', sprintf('%.3E', norm(g_n))])
        end
    end
end

