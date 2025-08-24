function Lv = Lv_evap(i,T,T_inf,fuel_names)
    gota  = clase_gota(fuel_names(i));
    T_ref = (2/3)*T + (1/3)*T_inf;
    Lv    = calcula_Lv_fuel(gota,[],T_ref,[])';
end