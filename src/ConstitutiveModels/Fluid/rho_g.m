function rho_g = rho_g(ii,T,T_inf,y,y_inf,model_g,fuel_names)
    y_cell       = cell(model_g.nInerts+1,1);
    for i=1:model_g.nInerts+1
        y_cell{i} = y(i);
    end 
    T_ref        = (2/3)*T + (1/3)*T_inf;
    y_ref        = cell(model_g.nInerts+1,1);
    for i=1:model_g.nInerts+1
        y_ref{i} = (2/3)*y_cell{i} + (1/3)*y_inf{i};
    end
    Yi_eval      = zeros(model_g.nInerts, size(T_ref,2)); 
    for i=1:model_g.nInerts
        Yi_eval(i,:) = y_ref{i}(:)';
    end
    Yf_eval      = zeros(1, size(T_ref,2));
    gota         = clase_gota(fuel_names(ii));
    matrix       = Polyfit_properties_pureCompounds(gota, model_g.comp_inerts, model_g.Tmin, model_g.Tmax, model_g.P, model_g.N_polyfit);
    Yf_eval(1,:) = 1-sum(Yi_eval);
    rho_g        = MixtureRules('rho_gas', T_ref, Yi_eval, Yf_eval, matrix, gota, model_g.comp_inerts, model_g.P);
end