function value = calcula_Pvap_fuel(obj, ~, Teval,~)

% Función que calcula la presión de vapor de cada componente puro de la
% mezcla. Análogamente a calcula_Lv no se calculan valores promedio, sino
% que interesa diferenciar entre componentes.

N_comb=length(obj.tipo_combustible);

for id_fuel=1:N_comb
    value(id_fuel)=devuelve_propiedad_fuel(obj, 'Pvap',obj.tipo_combustible(id_fuel),'liquido',Teval);
end

end
