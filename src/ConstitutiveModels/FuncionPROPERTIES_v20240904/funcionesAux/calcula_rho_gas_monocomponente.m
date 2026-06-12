function rho = calcula_rho_gas_monocomponente(tipo_gas, T, P, R, gota)
%%Devuelve rho del gas especificado a la temperatura T y presion P
% GAS IDEAL. Abierto a añadir gas real en futura ampliación.

if strcmpi(tipo_gas,'N2')
   mw = 28/1000 ; %(* kg/mol *) 
elseif strcmpi(tipo_gas,'O2')
   mw = 32/1000 ; %(* kg/mol *)
elseif strcmpi(tipo_gas,'CO2')
   mw = 44/1000;  %(* kg/mol *)    
elseif strcmpi(tipo_gas,'H2O')
   mw = 18/1000 ;  %(* kg/mol *) 
elseif strcmpi(tipo_gas,'CO')
   mw = 28/1000 ;  %(* kg/mol *)
else %(it is a fuel)
   index = find(strcmp(gota.tipo_combustible, tipo_gas));
   mw=gota.mw(index);
end

rho=P*mw/(R*T);


end

