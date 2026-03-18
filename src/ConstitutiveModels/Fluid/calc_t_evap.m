function t_evap = calc_t_evap(T_inf,y_inf,model_l,model_g,frac_masL,fuel_names,x_lg)
    % Calculate T_i
    T_i = calc_T_i(T_inf,y_inf,model_l,model_g,fuel_names,frac_masL);
    
    % Calculate t_evap
    if model_l.nSpecies==1
        T_ref       = (2/3)*T_i + (1/3)*T_inf;
        y_L         = frac_masL';
        y_R         = calc_y_gInter(y_L,model_l,model_g,T_i);
        y_ref       = cell(model_g.nSpecies,1);
        for ii=1:model_g.nSpecies
            y_ref{ii} = (2/3)*y_R{ii} + (1/3)*y_inf{ii};
        end

        polynomial = squeeze(model_l.matrix{1,model_g.nInerts+1});
        rho_l_evap = Polyval(polynomial, T_i);
        rho_g_evap = calc_rho(model_g,y_ref,T_ref);
        D_T_evap   = calc_D_T(T_ref,y_ref,model_g);
        Cp_evap    = calc_Cp(T_ref,y_ref,model_g);
        Lv         = calcula_Lv_fuel(model_l.gota,[],T_ref,[]);
        t_evap     = ((2*x_lg)^2 * rho_l_evap) / (D_T_evap * rho_g_evap * 8 * log(1 + (Cp_evap * (T_inf - T_i)) / Lv));
        t_evap     = real(t_evap);
    elseif model_l.nSpecies>1
        t_evap_i   = calc_t_evap_i(T_i,T_inf,y_inf,model_l,model_g,frac_masL,fuel_names,x_lg);
        t_evap     = sum(t_evap_i);
        t_evap     = real(t_evap);
    else
        disp('Error calculating t_evap')
        return
    end
end