function Pv = calc_Pvap(model_l,T)
    fuel = model_l.gota;
    Pv   = zeros(model_l.nSpecies,1);
    for II=1:model_l.nSpecies
        Pv(II) = devuelve_propiedad_fuel(fuel, 'pvap', fuel.tipo_combustible{II}, 'liquido', T);
    end
end