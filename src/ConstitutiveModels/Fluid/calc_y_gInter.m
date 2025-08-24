function y_gInter = calc_y_gInter(y_l, model_l, model_g, T_g)

    y_gInter    = cell(model_g.nSpecies,1);
    
    for i=1:model_l.nSpecies
        y_l{i}      = y_l{i}(end,end);
    end
    
    fuel        = model_l.gota;
    MW_g        = [model_g.Inerts_mw fuel.mw]';
    
    Pv          = calc_Pvap(model_l,T_g(1));
    
    sumY_W      = zeros(1,1);
    for i=1:model_l.nSpecies
        sumY_W = sumY_W + y_l{i}/MW_g(i+model_g.nInerts);
    end
    
    X_l_fuels   = zeros(model_l.nSpecies,1);
    P_i         = zeros(model_l.nSpecies,1);
    for i=1:model_l.nSpecies
        X_l_fuels(i) = (y_l{i}/MW_g(i+model_g.nInerts))/sumY_W;
        P_i(i)       = X_l_fuels(i)*Pv(i);
    end
    
    X_g_fuels   = P_i./model_g.P;
    
    X_g_Tinerts = 1-sum(X_g_fuels);
    X_g_Inerts  = zeros(model_g.nInerts,1);
    for i=1:model_g.nInerts
        X_g_Inerts(i)  = X_g_Tinerts*model_g.frac_masG{i};
    end
    
    X_g         = [X_g_Inerts; X_g_fuels];
    
    sumXW       = zeros(1,1);
    for i=1:model_g.nSpecies
       sumXW = sumXW + X_g(i)*MW_g(i);
    end
    
    for i=1:model_g.nSpecies
        y_gInter{i} = (X_g(i)*MW_g(i))/sumXW;
    end

end