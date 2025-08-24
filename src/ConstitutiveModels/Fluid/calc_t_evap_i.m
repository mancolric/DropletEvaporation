function t_evap_i = calc_t_evap_i(T_i,T_inf,y_inf,model_l,model_g,frac_masL,fuel_names,x_lg)
    t_evap_i = zeros(model_l.nSpecies,1);
    for ii=1:model_l.nSpecies    
        T_ref      = (2/3)*T_i(ii) + (1/3)*T_inf;
        y_L        = frac_masL';
        y_R        = calc_y_gInter(y_L,model_l,model_g,T_i(ii));
        y_ref      = cell(model_g.nSpecies,1);
        for i=1:model_g.nSpecies
            y_ref{i} = (2/3)*y_R{i} + (1/3)*y_inf{i};
        end
        
        Yi_eval    = zeros(model_g.nInerts, size(T_ref,2)); 
        for II=1:model_g.nInerts
            Yi_eval(II,:) = y_ref{II}(:)';
        end
        Yf_eval    = zeros(1, size(T_ref,2));
        V          = zeros(model_l.nSpecies,1);
        V(1)       = 4/3*pi*x_lg^3;
        V(end)     = V(1)*frac_masL{end};
        a          = zeros(model_l.nSpecies,1);
        a(end)     = (V(end)*3/(4*pi))^(1/3);
        a(1)       = x_lg; 
        polynomial   = squeeze(model_l.matrix{1,model_g.nInerts+ii});
        rho_l_evap   = Polyval(polynomial, T_i(ii));
        polynomial   = squeeze(model_g.matrix{6,model_g.nInerts+ii});
        Cp_gas       = Polyval(polynomial, T_ref);
        
        gota         = clase_gota(fuel_names(ii));
        matrix       = Polyfit_properties_pureCompounds(gota, model_g.comp_inerts, model_g.Tmin, model_g.Tmax, model_g.P, model_g.N_polyfit);
        Yf_eval(1,:) = 1-sum(Yi_eval);
        rho_g_evap   = MixtureRules('rho_gas', T_ref, Yi_eval, Yf_eval, matrix, gota, model_g.comp_inerts, model_g.P);
        k_g_evap     = MixtureRules('k_gas', T_ref, Yi_eval, Yf_eval, matrix, gota, model_g.comp_inerts, model_g.P);
        D_T_evap     = k_g_evap./(rho_g_evap.*Cp_gas(:)');
        Lv           = calcula_Lv_fuel(gota,[],T_ref,[]);
        
        da2_dt       = 2*(rho_g_evap/rho_l_evap)*D_T_evap*log(1+(Cp_gas*(T_inf-T_i(ii)))/Lv);
        if ii==model_l.nSpecies
            t_evap_i(ii) = a(ii)^2/da2_dt;
        else
            t_evap_i(ii) = (a(ii)^2-a(ii+1)^2)/da2_dt;
        end
    end
end