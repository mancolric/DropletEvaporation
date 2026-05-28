function ap = calc_ap(T_m, y, model)

[m,n]   = size(T_m);
T_v     = T_m(:)'; % T_v[1 x nElems*(2*p+1)]

% Compositions at the nodes: 2 matrices, one for inerts (Yi_eval) another for fuels (Yf_eval)
% Inert order: N2 ; O2 ; CO2 ; H2O (the one stated in 'comp_inerts')
% Fuel order: the one stated in 'fuel_names'
% Each column is a node
% Each row is a compoundYi_eval = zeros(size(model.comp_inerts,2), size(T_v,2)); % Liq (Yi_eval[nInerts x nElems*(2*p+1)])

if model.bool_liq==1 % Liquid
    warning('¿Radiación en líquido ?')
    ap = zeros(size(T_m));
elseif model.bool_liq==0 % Gas
    Yi_eval = zeros(model.nInerts, size(T_v,2)); 
    for II=1:model.nInerts
    Yi_eval(II,:) = y{II}(:)';
    end

    Yf_eval = zeros(model.nSpecies-model.nInerts, size(T_v,2)); 
    for II=1:model.nSpecies-model.nInerts

        Yf_eval(II,:) = y{II+model.nInerts}(:)';
    end
    cells_ap = MixtureRules('ap_gas', T_v, Yi_eval, Yf_eval, model.matrix, model.gota, model.comp_inerts, model.P);
    ap   = cell(length(cells_ap),1);
    for II=1:model.nSpecies

        ap{II} = reshape(cells_ap{II},m,n); % D_rho[nElems x (2*p+1)]

    end
else
    error('Error calc_ap')
end

end