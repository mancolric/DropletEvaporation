function T = calc_T(H_m, rhoy, model)

[rho_m,y] = calc_rho_y(rhoy,model);
[m,n]     = size(rho_m);

% Initial T
if model.bool_liq==1
    T0      = repmat(300, m*n, 1);
else
    T0      = repmat(350, m*n, 1);
end

%Residual:
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

%Preconditioned residual:
[~,J]       = Residual(T0, true);
Jv          = diag(J);
function g = PrecResidual(T)
    [f,~]   = Residual(T, false);
    g       = Jv.\f;
end

% [T, Iteration, flag]    = NewtonRaphsonLS(@Residual, T0, 0.0, 1e-8, 100, 100);
[T, Iteration, flag]    = NewtonRaphson(@Residual, T0, 0.0, 1e-8, 10);
% [T, Iteration, flag]    = Anderson(@PrecResidual, T0, 1e-8, 100, 0);

T = reshape(T,m,n); % T[nElems x (2*p+1)] = T_m

if true && flag<0
    warning('Evaluation of temperature failed')
    T               = T*NaN;    
elseif flag<0
%     error("Error en el cálculo de T, Máximo iteraciones Newton-Raphson")
    Tv              = T(:);
    [f,~]           = Residual(Tv, false);
    [fmax,imax]     = max(abs(f));
    display(Tv(imax))
    
    %Plot h vs T:
    Tv              = linspace(10.0, 4000, 1000);
    y2              = cell(size(y));
    rhoy2           = cell(size(y));
    for ii=1:length(y)
        y2{ii}      = y{ii}(imax);
        rhoy2{ii}   = rhoy{ii}(imax);
    end
    display(y2)
    display(rhoy2)
    hv              = calc_h(Tv, y2, model);
    figure()
    plot(Tv, rho_m(imax)*hv-H_m(imax))
    
    T       = T*NaN;
end

end

