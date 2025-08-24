% ##########################################################################
% Código (ejemplo) donde se muestra cómo manejar propiedades para el casoR
% de mezclas, tanto en fase líquida como en fase gas. Este script deberá
% ser integrado en el modelo de evaporación de gota de la forma que consideréis más óptima.
% ##########################################################################

clear
clc

% addpath(genpath(pwd))


% ###########################################################################################################
% Step 1: Adjust individual properties of all the compounds (fuels +
% inerts) to a polynomial of degree N_polyfit. 
% NOTE: This only needs to be done ONCE for a given simulation.

% Name(s) of the fuel(s) contained in the  droplet:
fuel_names={'Metanol' 'Butanol'};
% fuel_names={'Etanol'};
gota=clase_gota(fuel_names);  % Creates object with constant properties: Tb, Tc, sigma, e_k, etc.


% Name(s) of the inert gases:
comp_inerts={'N2' 'O2' 'CO2' 'H2O'}; 

% Min. and max. temperatures to be evaluated in the simulation (K):
Tmin=300; Tmax=1700;

% Pressure (Pa):
P=101325;

% Build a matrix with the adjusted Polyfit coefficients for all the properties:
N_polyfit=5;
MatrixProperties=Polyfit_properties_pureCompounds(gota, comp_inerts, Tmin, Tmax, P, N_polyfit); % N_properties x N_inerts+N_fuels x poly_order+1

% ###########################################################################################################

% Step 2: apply mixture rules for obtaining the properties for any given
% mixture and temperature. This will be done many times during a droplet
% evaporation simulation.

% Temperatures at the node(s):
T_eval=linspace(300, 1700, 10); % T at node 1, 2, ..., N (Gas)
% T_eval=linspace(300, 350, 10); % T at node 1, 2, ..., N (Liq)

% Compositions at the nodes: 2 matrices, one for inerts (Yi_eval) another
% for fuels (Yf_eval)
% Inert order: N2 ; O2 ; CO2 ; H2O (the one stated in 'comp_inerts')
% each column is a node
% each row is a compound
% [Example: Gas containing 2 inerts: N2 and CO2. 10 nodes. Decreasing mass fractsions. For any given node, Y_N2/Y_CO2=3]
Yi_eval= [0.75*linspace(0.01,0.99, 10); 0*linspace(0.01,0.99, 10); 0.25*linspace(0.01,0.99, 10); 0*linspace(0.01,0.99, 10)]; % Y at node 1, 2, ..., N for all inerts


% Fuel order: the one stated in 'fuel_names'
% each column is a node
% each row is a compound
% [Example: Gas containing 2 fuels (the ones declared in fuel_names). Same 10 nodes as for the inerts] 
Yf_eval(1,:)=(ones(1, size(Yi_eval,2))-sum(Yi_eval))*0.6;
Yf_eval(2,:)=(ones(1, size(Yi_eval,2))-sum(Yi_eval))*0.4;
% [CHECK for the gas]: sum(Yi_eval,1) + sum(Yf_eval,1) should equal 1 for
% all nodes


% Functions calculating any given property in terms of T_eval and Y_eval.
% T_eval is an array with N temperatures (for N nodes).
% Y_eval has two matrices: Yi_eval (inert gases), Yf_eval (fuels). Each matrix
% has N rows (N nodes) and M columns (M compounds).
liquid_property=0;
if liquid_property==1

    if size(fuel_names,2)==1 %Monocomponent
    % Liquid properties for monocomponent case: no dependence with
    % composition, only with T. No mixture rules needed.
    row_rho_liq=1; row_C_liq=2; row_k_liq=3; id_f=5;
    vector_rho_liq = Polyval(squeeze(MatrixProperties{row_rho_liq, id_f}), T_eval);
    vector_C_liq = Polyval(squeeze(MatrixProperties{row_C_liq, id_f}), T_eval);
    vector_k_liq = Polyval(squeeze(MatrixProperties{row_k_liq, id_f}), T_eval);
    vector_D_liq = 0; % Mass diffusivity not required for a pure compound!
    vector_Dt_liq = vector_k_liq./(vector_rho_liq.*vector_C_liq); % Liquid thermal diffusivity Dt=k/(rho*C)

    else %Multicomponent
    % Liquid properties: no inerts inside the droplet --> Yi_eval substituted
    % by matrix of zeros (avoids potential input-related errors)
    Yi_liquid = zeros(size(comp_inerts,2), size(T_eval,2));
    vector_rho_liq = MixtureRules('rho_liq', T_eval, Yi_liquid, Yf_eval, MatrixProperties, gota, comp_inerts, P);
    vector_C_liq = MixtureRules('C_liq', T_eval, Yi_liquid, Yf_eval, MatrixProperties, gota, comp_inerts, P);
    vector_k_liq = MixtureRules('k_liq', T_eval, Yi_liquid, Yf_eval, MatrixProperties, gota, comp_inerts, P);
    vector_D_liq = MixtureRules('D_liq', T_eval, Yi_liquid, Yf_eval, MatrixProperties, gota, comp_inerts, P);
    vector_Dt_liq = vector_k_liq./(vector_rho_liq.*vector_C_liq); % Liquid thermal diffusivity Dt=k/(rho*C)
    end

else % Gas properties: both Yf_eval and Yi_eval required:
    % The gas is always a mixture, even if the droplet only contains 1
    % compound. Therefore we always use mixture rules.
    vector_rho_gas = MixtureRules('rho_gas', T_eval, Yi_eval, Yf_eval, MatrixProperties, gota, comp_inerts, P);
    vector_Cp_gas = MixtureRules('Cp_gas', T_eval, Yi_eval, Yf_eval, MatrixProperties, gota, comp_inerts, P);
    vector_k_gas = MixtureRules('k_gas', T_eval, Yi_eval, Yf_eval, MatrixProperties, gota, comp_inerts, P);
    vector_D_gas = MixtureRules('D_gas', T_eval, Yi_eval, Yf_eval, MatrixProperties, gota, comp_inerts, P);
    vector_Dt_gas = vector_k_gas./(vector_rho_gas.*vector_Cp_gas); % Gas thermal diffusivity Dt=k/(rho*C)

end

