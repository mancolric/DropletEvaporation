%[x_n, Iteration, flag] = Anderson(fun, x0, TolG, MaxIter, m)
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
%accelerated with Anderson's method. "m" is the number of last iterations
%saved.
%
%fun receives the argument x and returns g(x)
function [x_n, nIters, flag] = Anderson(fun, x0, TolG, MaxIter, m)

    %First iteration:
    nIters      = 0;
    x_n         = x0;                   %solution
    g_n         = fun(x_n);             %preconditioned residual    
    p_n         = -g_n;                 %step
    m           = min(length(x0), m);
    Xm          = zeros(length(x0),m);  %memory matrices
    Gm          = zeros(length(x0),m);
    mXG         = 0;                    %nb of saved solutions
    while true

        %Exit loop:
        if norm(g_n)<=TolG
            flag    = 1;
            break
        elseif nIters==MaxIter
            warning('Reached maximum of iterations')
            flag   = -1;
            break
        elseif any(isnan(g_n))
            warning('Found NaNs in preconditioned residual')
            flag    = -2;
            break
        end

        %Solution and preconditioned residual at next iteration:
        x_np1       = x_n + p_n;
        g_np1       = fun(x_np1);
        
        %Update X and G matrices with variations in X and G:
        mXG         = min(m, mXG+1);
        Deltag_norm = norm(g_np1-g_n) + 1e-12;
        Xm(:,2:mXG) = Xm(:,1:mXG-1);
        Xm(:,1)     = (x_np1 - x_n)/Deltag_norm;
        Gm(:,2:mXG) = Gm(:,1:mXG-1);
        Gm(:,1)     = (g_np1 - g_n)/Deltag_norm;
        
        %Update (x_n, g_n):
        nIters      = nIters+1;
        x_n         = x_np1;
        g_n         = g_np1;
        
        %#Next step. Apply multisecant formula:
        %   pn  = - H g^n = - [I + (X-G) (G^T G)^{-1} G^T] g^n 
        %       = -g^n - X*gamma + G*gamma,
        %with gamma = (G^T G)^{-1} G^T g^n 
        %If G=U Sigma V^T, then
        %   gamma = V Sigma^{-1} U^T g^n
        %Note that Anderson acceleration would read
        %   x_star      = xn - X*gamma
        %   g_star      = gn - G*gamma
        %   x_(n+1)     = x_star - g_star = xn - X*gamma - (gn - G*gamma)
        %hence
        %   p_n         = - gn - X*gamma + G*gamma
        if any(isnan(g_np1)) || ~all(isfinite(g_np1))
            warning('NaNs or Infs in preconditioned residual')
            flag        = -2;
            break
        end
        [U,S,V]                 = svd(Gm,0);
        S_diag                  = diag(S);
        S_diag(S_diag<1e-10)    = Inf;
        S_inv                   = diag(1.0./S_diag); 
        gamma                   = V*(S_inv*(transpose(U)*g_n));
        p_n                     = -g_n - Xm*gamma + Gm*gamma;
        
    end
    
end

