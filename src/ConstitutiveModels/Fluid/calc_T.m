function T = calc_T(H_m, rhoy, model)

[rho_m,y] = calc_rho_y(rhoy,model);
[m,n]     = size(rho_m);

% Initial T
if model.bool_liq==1
    T0      = repmat(300, m*n, 1);
else
    T0      = repmat(350, m*n, 1);
end

% NewtonRaphson: T from h, rho, y, H y Cp:
    % rho = sum(rho * y)
    % rho * h (T,y1,...,yN) = H
    % sum(rho * y) * h (T,y1,...,yN) = H
    % f (T) = 0
    function [f,J] = Residual(T_v, ComputeJ)
        T_m     = reshape(T_v,m,n);
        h       = calc_h(T_m,y,model);
        f       = rho_m.*h-H_m;
        f       = f(:);
        if ComputeJ
            Cp  = calc_Cp(T_m,y,model);
            J   = rho_m(:).*Cp(:);
            J   = sparse(1:length(f), 1:length(f), J);
        else
            J  = NaN;
        end
    end

% [T, Iteration, flag]    = NewtonRaphsonLS(@Residual, T0, 0.0, 1e-8, 100, 100);
[T, Iteration, flag]    = NewtonRaphson(@Residual, T0, 0.0, 1e-8, 100);

T = reshape(T,m,n); % T[nElems x (2*p+1)] = T_m

if flag==-2
    error("Error en el cálculo de T Máximo iteraciones LS")
elseif flag==-1
    error("Error en el cálculo de T, Máximo iteraciones Newton-Raphson")
end

end

