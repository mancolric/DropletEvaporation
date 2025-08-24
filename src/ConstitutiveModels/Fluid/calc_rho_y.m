function [rho_m,y] = calc_rho_y(rhoy,model)

rho_m = zeros(size(rhoy{1}));

for II=1:model.nSpecies
    rho_m = rho_m + rhoy{II}; % rho_m[nElems x (2*p+1)]
end

y = cell(model.nSpecies,1);
for II=1:model.nSpecies
    y{II} = rhoy{II}./rho_m; % y[nElems x (2*p+1)]
end
