function value = calcula_ReactionOrders(obj, bool_liq)

if bool_liq

    error("Arrhenius in liquid?")

else
N_comb=length(obj.fuel);
N_iner=length(obj.comp_inerts);

for id_gas=1:obj.nReaction
    value(id_gas,:)=devuelve_propiedad_fuel(obj.gota, 'ROrder',obj.species(N_comb+N_iner),'vapor',298.15);
end
end
end