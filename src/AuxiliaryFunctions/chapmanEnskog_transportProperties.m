function devuelve= chapmanEnskog_transportProperties(compuesto,gota,T,model)
% https://books.google.es/books?id=BVSH2AkjpP8C&pg=PA888&lpg=PA888&dq=chapman+enskog+thermal+conductivity+octane&source=bl&ots=VyzcHvF-w2&sig=ACfU3U124HrwU__VtUohizDWRM1hsLvf3Q&hl=es&sa=X&ved=2ahUKEwiNjLOt0J3iAhWlDGMBHZ71DA0Q6AEwCHoECAYQAQ#v=onepage&q&f=false
% load('Matriz_Omegas.mat')
% Funcion que calcula, para un compuesto puro (octano/aire) sus propiedades
% de transporte de acuerdo a la teoria cinetica de Chapman-Enskog
% (utilizada en paper Hubbard)

% OJO UNIDADES, ESPECIALMENTE MW!
if strcmp(compuesto,'octano')
% PARA OCTANO:
sigma=7.407;
epsilon_k=333;
MW=gota.mw*1000; % g/mol
Cpv=devuelve_propiedad_fuel(gota, 'cp', 'Octano', 'vapor', T ); %( %(* J/kg/K *)Si se evalúa propiedad de transporte gas, gota ya vendrá actualizada con T_1/3)

elseif strcmp(compuesto,'aire')
% PARA AIRE:
sigma=3.617;
epsilon_k=97;
MW=28.97;  % g/mol
Cpv=0.233*calcula_Cp_gas_monocomponente('O2', T)+0.767*calcula_Cp_gas_monocomponente('N2', T); %Cp de aire  %(* J/kg/K *)

end


kT_e=T/epsilon_k; % Parametro para entrar a la tabla de las omegas
omega_mu=model.F1(kT_e); % DEvuelve omega_mu

% Viscosidad:
mu_v=2.67e-6*((MW*T)^0.5)/(sigma^2*omega_mu); %(m^2/s)

% Conductividad:
omega_lambda=omega_mu;
lambda_v=8.32e-2*((T/MW)^0.5)/(sigma^2*omega_lambda) + 1.32*mu_v*(Cpv - 5/2*8.314/(MW/1000)); %Ojo! En segundo termino MW debe ir en kg/mol (por unidades cte gases ideales)


% Coeficiente D_12 (igual para ambos):
MW_1=gota.mw*1000; % g/mol;
MW_2=28.97; % g/mol;
sigma_12=0.5*(3.617+7.407);
epsilon_k=(97*333)^0.5;
kT_e=T/epsilon_k; % Parametro para entrar a la tabla de las omegas
omega_D=model.F2(kT_e); % DEvuelve omega_D
Presion=1; %atm

D_12=1.86e-7*((T^3*(1/MW_1+1/MW_2))^0.5)/(sigma_12^2*omega_D*Presion);


devuelve=[mu_v,lambda_v, D_12];
end
