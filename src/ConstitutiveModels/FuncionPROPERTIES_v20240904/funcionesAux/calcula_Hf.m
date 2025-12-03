function value = calcula_Hf(obj, bool_liq)

% Realmente NO se calcula el Lv de una mezcla, sino que se evalúan por
% separado cada uno de los valores de Lv de cada componente. Se utiliza
% esta función como intermediaria entre actualiza_gota y devuelve_propiedad
% tan solo por mantener la coherencia con respecto al resto de
% propiedades:

if bool_liq
N_comb=length(obj.tipo_combustible);

for id_fuel=1:N_comb
    value(id_fuel)=devuelve_propiedad_fuel(obj, 'Hf',obj.tipo_combustible(id_fuel),'vapor',298.15);
end

else
N_comb=length(obj.fuel);
N_iner=length(obj.comp_inerts);

for id_gas=1:N_iner
    value(id_gas)=devuelve_propiedad_gas(obj, 'Hf',obj.species(id_gas),'vapor',298.15);
end

for id_gas=N_iner+1:N_iner+N_comb
    value(id_gas)=devuelve_propiedad_fuel(obj.gota, 'Hf',obj.species(id_gas),'vapor',298.15);
end

end
end