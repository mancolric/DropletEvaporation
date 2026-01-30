function value = calcula_fuel_mw(obj, bool_liq)

if bool_liq
N_comb=length(obj.tipo_combustible);

for id_fuel=1:N_comb
    value(id_fuel)=devuelve_propiedad_fuel(obj, 'MW',obj.tipo_combustible(id_fuel),'vapor',298.15);
end

else
N_comb=length(obj.fuel);
N_iner=length(obj.comp_inerts);

for id_gas=N_iner+1:N_iner+N_comb
    value(id_gas-N_iner)=devuelve_propiedad_fuel(obj.gota, 'MW',obj.species(id_gas),'vapor',298.15);
end

end
end