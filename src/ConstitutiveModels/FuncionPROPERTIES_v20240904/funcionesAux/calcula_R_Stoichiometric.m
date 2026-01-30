function value = calcula_R_Stoichiometric(obj, bool_liq)

if bool_liq

    error("Arrhenius in liquid?")

else
N_comb=length(obj.fuel);
N_iner=length(obj.comp_inerts);

for id_gas=N_iner+1:N_iner+N_comb
    value(id_gas,:)=devuelve_propiedad_fuel(obj.gota, 'Rsoi',obj.species(id_gas),'vapor',298.15);
end
end
end