function value = calcula_k_fuel(obj, ~, Teval,~)

% Realmente NO se calcula el Lv de una mezcla, sino que se evalúan por
% separado cada uno de los valores de Lv de cada componente. Se utiliza
% esta función como intermediaria entre actualiza_gota y devuelve_propiedad
% tan solo por mantener la coherencia con respecto al resto de
% propiedades:

N_comb=length(obj.tipo_combustible);

for id_fuel=1:N_comb
    value(id_fuel)=devuelve_propiedad_fuel(obj, 'lambda',obj.tipo_combustible(id_fuel),'liquido',Teval);
end

end