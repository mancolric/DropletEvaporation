function value = calcula_Arrhenius(obj, bool_liq)

if bool_liq

    error("Arrhenius in liquid?")

else
N_comb=length(obj.fuel);
N_iner=length(obj.comp_inerts);

for id_gas=1:N_iner
    value(id_gas,:)=[0.0, 0.0, 0.0];
end

for id_gas=N_iner+1:N_iner+N_comb
    value(id_gas,:)=devuelve_propiedad_fuel(obj.gota, 'Arr',obj.species(id_gas),'vapor',298.15);
end

end
end