function h_i = calc_h_i(T_m,model)

h_i   = cell(model.nSpecies,1);
T_ref = 298.15;
for i=1:model.nSpecies
    if model.bool_liq==1
        coeff     = model.matrix{2,i+model.nInerts}.p;
        mu        = model.matrix{2,i+model.nInerts}.mu; % mu(1)=mean | mu(2)=std
        N         = length(coeff);
        div       = N:-1:1;
        coef_p    = [coeff(:)' ./ div, 0];

        T_n       = (T_m-mu(1))/mu(2);
        T_ref_n   = (T_ref-mu(1))/mu(2);

        % h_i_T     = (T_n.*(coeff(6) + T_n.*(coeff(5)/2 + T_n.*(coeff(4)/3 + T_n.*(coeff(3)/4 + T_n.*(coeff(2)/5 + T_n.*(coeff(1)/6)))))));
        % h_i_Tref  = (T_ref_n.*(coeff(6) + T_ref_n.*(coeff(5)/2 + T_ref_n.*(coeff(4)/3 + T_ref_n.*(coeff(3)/4 + T_ref_n.*(coeff(2)/5 + T_ref_n.*(coeff(1)/6)))))));
        h_i_T     = polyval(coef_p, T_n);
        h_i_Tref  = polyval(coef_p, T_ref_n);

        h_i{i}    = (h_i_T-h_i_Tref).*mu(2) - model.Lv(:,i) + model.Hf(:,i);
    else
        coeff     = model.matrix{6,i}.p;
        mu        = model.matrix{6,i}.mu; % mu(1)=mean | mu(2)=std
        N         = length(coeff);
        div       = N:-1:1;
        coef_p    = [coeff(:)' ./ div, 0];

        T_n       = (T_m-mu(1))/mu(2);
        T_ref_n   = (T_ref-mu(1))/mu(2);

        % h_i_T     = (T_n.*(coeff(6) + T_n.*(coeff(5)/2 + T_n.*(coeff(4)/3 + T_n.*(coeff(3)/4 + T_n.*(coeff(2)/5 + T_n.*(coeff(1)/6)))))));
        % h_i_Tref  = (T_ref_n.*(coeff(6) + T_ref_n.*(coeff(5)/2 + T_ref_n.*(coeff(4)/3 + T_ref_n.*(coeff(3)/4 + T_ref_n.*(coeff(2)/5 + T_ref_n.*(coeff(1)/6)))))));
        h_i_T     = polyval(coef_p, T_n);
        h_i_Tref  = polyval(coef_p, T_ref_n);

        h_i{i}    = (h_i_T-h_i_Tref).*mu(2) + model.Hf(:,i);
    end
end

end