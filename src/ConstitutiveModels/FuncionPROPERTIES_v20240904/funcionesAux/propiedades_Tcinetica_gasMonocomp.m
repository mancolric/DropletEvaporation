function vector_propiedades= propiedades_Tcinetica_gasMonocomp(tipo_gas)
% Función que me devuelve las propiedades de teoría cinética para un gas
% monocomponente dado

if strcmpi(tipo_gas,'N2')
   mw = 28/1000 ; %(* kg/mol *)
   sigma=3.621*1e-10; % (m) Lennard-Jones length N2, CHEMKIN transport
   e_k=97.53; % (K) Lennard-Jones energy N2, CHEMKIN transport
   Zrot_298=4; % Sacado de pdf transporte de CHEMKIN
   molecula_lineal=0; 
   
elseif strcmpi(tipo_gas,'O2')
   mw = 32/1000 ; %(* kg/mol *)
   sigma=3.458*1e-10; % (m) Lennard-Jones length N2, CHEMKIN transport
   e_k=107.4; % (K) Lennard-Jones energy N2, CHEMKIN transport
   Zrot_298=3.8; % Sacado de pdf transporte de CHEMKIN
   molecula_lineal=0; % Sacado de pdf transporte de CHEMKIN
    
elseif strcmpi(tipo_gas,'CO2')
   mw = 44/1000;  %(* kg/mol *)
   sigma=3.763*1e-10; % (m) Lennard-Jones length N2, CHEMKIN transport
   e_k=244; % (K) Lennard-Jones energy N2, CHEMKIN transport
   Zrot_298=2.1; % Sacado de pdf transporte de CHEMKIN
   molecula_lineal=0; % Sacado de pdf transporte de CHEMKIN
        
elseif strcmpi(tipo_gas,'H2O')
   mw = 18/1000 ;  %(* kg/mol *) 
   sigma=2.605*1e-10; % (m) Lennard-Jones length N2, CHEMKIN transport
   e_k=572.4; % (K) Lennard-Jones energy N2, CHEMKIN transport
   Zrot_298=4; % Sacado de pdf transporte de CHEMKIN
   molecula_lineal=0; % Sacado de pdf transporte de CHEMKIN

elseif strcmpi(tipo_gas,'CO')
   mw = 28/1000 ; %(* kg/mol *)
   sigma=3.650*1e-10; % (m) Lennard-Jones length CO, CHEMKIN transport
   e_k=98.100; % (K) Lennard-Jones energy CO, CHEMKIN transport
   Zrot_298=1.800; % Sacado de pdf transporte de CHEMKIN
   molecula_lineal=1;

end


vector_propiedades=[mw, sigma, e_k, Zrot_298, molecula_lineal];

return

