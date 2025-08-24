function [ec,verify] = calc_ecPvap(y_l, y_g, model_l, model_g, T)

fuel                    = model_l.gota;
MW_l                    = fuel.mw';
MW_g                    = [model_g.Inerts_mw fuel.mw]';

sum_L                   = 0;
for II=1:model_l.nSpecies
    sum_L               = sum_L + y_l{II}/MW_l(II);
end
sum_G                   = 0;
for II=1:model_g.nSpecies
    sum_G               = sum_G + y_g{II}/MW_g(II);
end

Pv                      = calc_Pvap(model_l,T);

sumYW                   = 0;
for i=1:model_l.nSpecies
    sumYW = sumYW + y_l{i}/MW_g(i+model_g.nInerts);
end

X_l_fuels               = zeros(model_l.nSpecies,1);
P_i                     = zeros(model_l.nSpecies,1);
for i=1:model_l.nSpecies
    X_l_fuels(i) = (y_l{i}/MW_g(i+model_g.nInerts))/sumYW;
    P_i(i)       = X_l_fuels(i)*Pv(i);
end

ec                      = zeros(model_l.nSpecies,1);
for II=1:model_l.nSpecies
    ec(II) = (y_g{II+model_g.nInerts}/MW_g(II+model_g.nInerts))*sum_L - (Pv(II)/model_l.P)*(y_l{II}/MW_l(II))*sum_G;
end

X_g                     = zeros(model_l.nSpecies,1);
verify                  = zeros(model_l.nSpecies,1);
for II=1:model_l.nSpecies
    X_g(II)     = (y_g{II+model_g.nInerts}/MW_g(II+model_g.nInerts))/sum_G;
    verify(II)  = (P_i(II)-model_l.P*X_g(II))/P_i(II);
end

end