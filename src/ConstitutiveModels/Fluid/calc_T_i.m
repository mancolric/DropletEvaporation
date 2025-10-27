function T_i = calc_T_i(T_inf,y_inf,model_l,model_g,fuel_names,frac_masL)
    Nu    = 2;
    T_i   = zeros(model_l.nSpecies,1);
    for i=1:model_l.nSpecies
        fuel    = clase_gota(fuel_names(i));
        T_ini   = fix(0.8*fuel.Tb);
        eq_Ti   = @(Ti) Nu * lambda_g(i, Ti, T_inf, Y_vi(i, Ti, model_g, fuel_names, frac_masL), y_inf, model_g, fuel_names) * (T_inf - Ti) ...
                - (2 * rho_g(i, Ti, T_inf, Y_vi(i, Ti, model_g, fuel_names, frac_masL), y_inf, model_g, fuel_names)  ...
                * D_g(i, Ti, T_inf, Y_vi(i, Ti, model_g, fuel_names, frac_masL), y_inf, model_g, fuel_names) * Lv_evap(i, Ti, T_inf, fuel_names) ...
                * Y_vi_fuel(i, Ti, model_g, fuel_names, frac_masL) / (1 - Y_vi_fuel(i, Ti, model_g, fuel_names, frac_masL)));
        options = optimoptions('fsolve','Display','off','TolFun', 1e-6);
        T_i(i) = fsolve(eq_Ti,T_ini,options);
    end
end