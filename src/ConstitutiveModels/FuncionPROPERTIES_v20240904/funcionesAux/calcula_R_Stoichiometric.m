function value = calcula_R_Stoichiometric(obj, bool_liq)

if bool_liq

    error("Arrhenius in liquid?")

else
N_comb=length(obj.fuel);
N_iner=length(obj.comp_inerts);

value(1:obj.nReaction,:)=devuelve_propiedad_fuel(obj.gota, 'Rsoi',obj.species(N_iner+N_comb),'vapor',298.15);

end
end