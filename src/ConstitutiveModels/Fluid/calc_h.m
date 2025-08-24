function [h,h_i] = calc_h(T_m, y, model)

h   = zeros(size(T_m));
h_i = calc_h_i(T_m,model);
for i=1:model.nSpecies
    h   = h + h_i{i}.*y{i};
end

end