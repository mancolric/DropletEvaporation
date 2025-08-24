function Yv_fuel = Y_vi_fuel(i,Ti,model_g,fuel_names,frac_masL)
    Yv          = zeros(1+model_g.nInerts,1);

    y_l         = frac_masL{i};
    
    fuel        = clase_gota(fuel_names(i));
    MW_g        = [model_g.Inerts_mw fuel.mw]';
    
    Pv          = devuelve_propiedad_fuel(fuel, 'pvap', fuel.tipo_combustible, 'liquido', Ti);
    
    sumY_W      =  y_l/MW_g(1+model_g.nInerts);
    
    X_l_fuels   = (y_l/MW_g(1+model_g.nInerts))/sumY_W;
    P_i         = X_l_fuels*Pv;
    
    X_g_fuels   = P_i./model_g.P;
    
    X_g_Tinerts = 1-sum(X_g_fuels);
    X_g_Inerts  = zeros(model_g.nInerts,1);
    for i=1:model_g.nInerts
        X_g_Inerts(i)  = X_g_Tinerts*model_g.frac_masG{i};
    end
    
    X_g         = [X_g_Inerts; X_g_fuels];
    
    sumXW       = zeros(1,1);
    for i=1:size(X_g)
       sumXW = sumXW + X_g(i)*MW_g(i);
    end
    
    for i=1:size(X_g)
        Yv(i) = (X_g(i)*MW_g(i))/sumXW;
    end
    Yv_fuel = Yv(1+model_g.nInerts);
end