function [m] = calc_m_ini(x_lg, rho_g, D_rho, y_R, y_inf, model_l)
    m       = zeros(1,1);
    D_rho_i = cell(model_l.nSpecies,1);
    for i=1:model_l.nSpecies
        D_rho_i{i} = D_rho.*y_R{i+model_l.nInerts};
    end
    for i=1+model_l.nInerts:model_l.nInerts+model_l.nSpecies
        m_i = rho_g.*D_rho_i{i-model_l.nInerts}.*log(1+(y_R{i}-y_inf{i})./(1-y_R{i}))./x_lg;
        m   = m + m_i;
    end
end