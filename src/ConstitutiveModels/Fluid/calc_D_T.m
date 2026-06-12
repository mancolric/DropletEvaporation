function D_T = calc_D_T(T_m, y, model)

[m,n]   = size(T_m);
T_v     = T_m(:)'; % T_v[1 x nElems*(2*p+1)]

% Compositions at the nodes: 2 matrices, one for inerts (Yi_eval) another for fuels (Yf_eval)
% Inert order: N2 ; O2 ; CO2 ; H2O (the one stated in 'comp_inerts')
% Fuel order: the one stated in 'fuel_names'
% Each column is a node
% Each row is a compoundYi_eval = zeros(size(model.comp_inerts,2), size(T_v,2)); % Liq (Yi_eval[nInerts x nElems*(2*p+1)])

if model.bool_liq==1 % Liquid
    Yi_eval = zeros(size(model.comp_inerts,2), size(T_v,2));
    Yf_eval = zeros(model.nSpecies, size(T_v,2));
    for II=1:model.nSpecies
        Yf_eval(II,:) = (ones(1, size(Yi_eval,2))-sum(Yi_eval)).*y{II}(:)';
    end
    if model.nSpecies==1 % Monocomponent
        row_rho_liq=1; row_C_liq=2; row_k_liq=3; id_f=model.nInerts+1;
        vector_rho = Polyval(squeeze(model.matrix{row_rho_liq, id_f}), T_v);
        vector_C   = Polyval(squeeze(model.matrix{row_C_liq, id_f}), T_v);
        vector_k   = Polyval(squeeze(model.matrix{row_k_liq, id_f}), T_v);
        vector_Dt  = vector_k./(vector_rho.*vector_C); % Dt = k / (rho*C)
    else % Multicomponent
        Yi_liquid  = zeros(size(model.comp_inerts,2), size(T_v,2));
        vector_rho = MixtureRules('rho_liq', T_v, Yi_liquid, Yf_eval, model.matrix, model.gota, model.comp_inerts, model.P);
        vector_C   = MixtureRules('C_liq', T_v, Yi_liquid, Yf_eval, model.matrix, model.gota, model.comp_inerts, model.P);
        vector_k   = MixtureRules('k_liq', T_v, Yi_liquid, Yf_eval, model.matrix, model.gota, model.comp_inerts, model.P);
        vector_Dt  = vector_k./(vector_rho.*vector_C); % Dt = k / (rho*C)
    end
elseif model.bool_liq==0 % Gas
    Yi_eval = zeros(model.nInerts, size(T_v,2)); 
    for II=1:model.nInerts
    Yi_eval(II,:) = y{II}(:)';
    end
    Yf_eval = zeros(model.nSpecies-model.nInerts, size(T_v,2)); 
    for II=1:model.nSpecies-model.nInerts
        Yf_eval(II,:) = y{II+model.nInerts}(:)';
    end
    vector_rho = MixtureRules('rho_gas', T_v, Yi_eval, Yf_eval, model.matrix, model.gota, model.comp_inerts, model.P);
    vector_Cp  = MixtureRules('Cp_gas', T_v, Yi_eval, Yf_eval, model.matrix, model.gota, model.comp_inerts, model.P);
    vector_k   = MixtureRules('k_gas', T_v, Yi_eval, Yf_eval, model.matrix, model.gota, model.comp_inerts, model.P);
    vector_Dt  = vector_k./(vector_rho.*vector_Cp(:)'); % Gas thermal diffusivity Dt=k/(rho*C)
else
    error('Error calc_D_T')
end

D_T = reshape(vector_Dt,m,n); % D_rho[nElems x (2*p+1)]

end