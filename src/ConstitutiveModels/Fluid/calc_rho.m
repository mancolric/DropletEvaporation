function rho = calc_rho(model, y, T_m)

[m,n]   = size(T_m);
T_v     = T_m(:)'; % T_v[1 x nElems*(2*p+1)]

% Compositions at the nodes: 2 matrices, one for inerts (Yi_eval) another for fuels (Yf_eval)
% Inert order: N2 ; O2 ; CO2 ; H2O (the one stated in 'comp_inerts')
% Fuel order: the one stated in 'fuel_names'
% Each column is a node
% Each row is a compoundYi_eval = zeros(size(model.comp_inerts,2), size(T_v,2)); % Liq (Yi_eval[nInerts x nElems*(2*p+1)])

if model.bool_liq==1 % Liquid
%     rho     = 711+0*T_m;
%     return
    Yi_eval = zeros(size(model.comp_inerts,2), size(T_v,2));
    Yf_eval = zeros(model.nSpecies, size(T_v,2));
    for II=1:model.nSpecies
        Yf_eval(II,:) = (ones(1, size(Yi_eval,2))-sum(Yi_eval)).*y{II}(:)';
    end
    if model.nSpecies==1 % Monocomponent
        row_rho_liq=1; id_f=5;
        vector_rho = Polyval(squeeze(model.matrix{row_rho_liq, id_f}), T_v);
    else % Multicomponent
        Yi_liquid  = zeros(size(model.comp_inerts,2), size(T_v,2));
        vector_rho = MixtureRules('rho_liq', T_v, Yi_liquid, Yf_eval, model.matrix, model.gota, model.comp_inerts, model.P);
    end
elseif model.bool_liq==0 % Gas
%     rho     = 0.2+0*T_m;
%     return
    Yi_eval = zeros(model.nInerts, size(T_v,2)); 
    for II=1:model.nInerts
        Yi_eval(II,:) = y{II}(:)';
    end
    Yf_eval = zeros(model.nSpecies-model.nInerts, size(T_v,2)); 
    for II=1:model.nSpecies-model.nInerts
        Yf_eval(II,:) = y{II+model.nInerts}(:)';
    end
    vector_rho = MixtureRules('rho_gas', T_v, Yi_eval, Yf_eval, model.matrix, model.gota, model.comp_inerts, model.P);
else
    error('Error calc_rho')
end

rho = reshape(vector_rho,m,n); % D_rho[nElems x (2*p+1)]

end