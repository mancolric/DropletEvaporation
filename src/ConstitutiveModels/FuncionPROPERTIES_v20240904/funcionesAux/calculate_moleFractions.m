function [Xi_eval, Xf_eval] = calculate_moleFractions(comp_inerts,Yi_eval, gota, Yf_eval)

% Calculates the Xi_eval and Xf_eval from the input mass fractions:
N_inerts=size(comp_inerts,2);
N_fuels=size(gota.tipo_combustible,2); 
for id_i=1:N_inerts
    props_kinetic_theory(id_i,:) = propiedades_Tcinetica_gasMonocomp(comp_inerts(id_i));
    MW_i(id_i, 1) = props_kinetic_theory(id_i,1); %kg/mol
    Xi_eval_int(id_i, :) = (Yi_eval(id_i,:)/MW_i(id_i));
end
for id_f=1:N_fuels
    MW_f(id_f, 1) = gota.mw(id_f); %kg/mol
    Xf_eval_int(id_f, :) = (Yf_eval(id_f,:)/MW_f(id_f));
end 
Xi_eval = Xi_eval_int./(sum(Xi_eval_int,1)+sum(Xf_eval_int,1));
Xf_eval = Xf_eval_int./(sum(Xi_eval_int,1)+sum(Xf_eval_int,1));


end

