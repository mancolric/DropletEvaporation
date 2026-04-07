function matrixProperties = Polyfit_properties_pureCompounds(gota, comp_inerts, T_min, T_max, Pressure, poly_order)
% Function that generates a matrix with the polynominal fitting
% coefficients for all pure compounds:
% --> Liquid properties: evaluated between Tmin and Tb
% --> Gas properties: evaluated between Tmin and Tmax

% The Matrix structure is as follows:
% rows: rho_l, C_l, k_l, mu_l, rho_g, Cp_g, k_g (more properties can be readily added; Dif. coeff. are treated separately)
% columns: Inert_1, Inert_2, Inert_3 ... Fuel_1, Fuel_2, Fuel_3...


% ###############################################
% tic
verification_plots=0; %1-activates visualization of the polynomial fittings

% INITIALIZE FINAL MATRIX PROPERTIES:
N_inerts=size(comp_inerts,2);
N_fuels=size(gota.tipo_combustible,2);
N_properties=7;
matrixProperties = cell(N_properties, N_inerts+N_fuels);

% INITIALIZE TEMPERATURE ARRAYS
Npoints=100;
rho_liq = zeros(1, Npoints);
C_liq = zeros(1, Npoints);
k_liq = zeros(1, Npoints);
mu_liq = zeros(1, Npoints);
rho_gas = zeros(1, Npoints);
Cp_gas = zeros(1, Npoints);
k_gas = zeros(1, Npoints);


% POLYFITS FOR INERT COMPOUNDS (not soluble in liquid --> only gas properties are required):
T_gas = linspace(T_min, T_max, Npoints);

for id_compound=1:N_inerts
    for id_T=1:size(T_gas,2)
    T_eval=T_gas(id_T);

    % Row 1: rho_l [inert]
    % Liq. property not required for inert compounds (array of 0)
    
    % Row 2: C_l [inert]
    % Liq. property not required for inert compounds (array of 0)

    % Row 3: k_l [inert]
    % Liq. property not required for inert compounds (array of 0)

    % Row 4: mu_l [inert]
    % Liq. property not required for inert compounds (array of 0)

    % Row 5: rho_g [inert]
    Ru=8.314; % J/mol/K
    % rho_gas(id_T) = calcula_rho_gas_monocomponente(comp_inerts{id_compound}, T_eval, Pressure, Ru, gota);
    rho_gas(id_T) = 1.0;

    % Row 6: Cp_g [inert]
    % Cp_gas(id_T) = calcula_Cp_gas_monocomponente(comp_inerts{id_compound}, T_eval);
    Cp_gas(id_T) = 1500;

    % Row 7: k_g [inert]
    % k_gas(id_T) = calcula_lambda_vapor_Kee(comp_inerts{id_compound}, 0, T_eval);
    k_gas(id_T) = 0.05;

   
    end
    % Fit the data to the polynomial and store the coefficients in the matrix:
    matrixProperties{5, id_compound}=Polyfit(T_gas, rho_gas, poly_order);
    matrixProperties{6, id_compound}=Polyfit(T_gas, Cp_gas, poly_order);
    matrixProperties{7, id_compound}=Polyfit(T_gas, k_gas, poly_order);


%     % Optional: Check Polyfit accuracy (disabled by default)
    if verification_plots==1
        figure(5)
        plot(T_gas, rho_gas, 'Marker', 'o')
        hold on
        plot(T_gas, Polyval(Polyfit(T_gas, rho_gas, poly_order), T_gas))
        title(['rho gas adjustement, order: ' num2str(poly_order)])

        figure(6)
        plot(T_gas, Cp_gas, 'Marker', 'o')
        hold on
        plot(T_gas, Polyval(Polyfit(T_gas, Cp_gas, poly_order), T_gas))
        title(['Cp gas adjustement, order: ' num2str(poly_order)])

        figure(7)
        plot(T_gas, k_gas, 'Marker', 'o')
        hold on
        plot(T_gas, Polyval(Polyfit(T_gas, k_gas, poly_order), T_gas))
        title(['k gas adjustement, order: ' num2str(poly_order)])

        
    end
end



% POLYFITS FOR FUELS (gas + liquid properties required):
for id_compound=N_inerts+1:N_inerts+N_fuels
% (T_liq depends on Tb; therefore, a different T_liq for each fuel)
% For each compound: the Tb varies. Update it:
T_liq = linspace(T_min, 510, Npoints);

    for id_T=1:size(T_gas,2)
    Tg_eval=T_gas(id_T);
    Tl_eval=T_liq(id_T);

    % Row 1: rho_l [fuel]
    rho_liq(id_T)=devuelve_propiedad_fuel(gota, 'rho', gota.tipo_combustible{id_compound-N_inerts}, 'liquido', Tl_eval);

    % Row 2: C_l [fuel]
    C_liq(id_T)=devuelve_propiedad_fuel(gota, 'cp', gota.tipo_combustible{id_compound-N_inerts}, 'liquido', Tl_eval);

    % Row 3: k_l [fuel]
    k_liq(id_T)=devuelve_propiedad_fuel(gota, 'lambda', gota.tipo_combustible{id_compound-N_inerts}, 'liquido', Tl_eval);
 
    % Row 4: mu_l [fuel]
    mu_liq(id_T)=devuelve_propiedad_fuel(gota, 'mu', gota.tipo_combustible{id_compound-N_inerts}, 'liquido', Tl_eval);

    % Row 5: rho_g [fuel]
    % rho_gas(id_T) = calcula_rho_gas_monocomponente(gota.tipo_combustible{id_compound-N_inerts}, Tg_eval, Pressure, Ru, gota);
    rho_gas(id_T) = 1.0;
    
    % Row 6: Cp_g [fuel]
    % Cp_gas(id_T) =devuelve_propiedad_fuel(gota, 'cp', gota.tipo_combustible{id_compound-N_inerts}, 'vapor', Tg_eval);
    Cp_gas(id_T) = 1500;
    
    % Row 7: k_g [fuel]
    % k_gas(id_T) = calcula_lambda_vapor_Kee(gota.tipo_combustible{id_compound-N_inerts}, gota, Tg_eval);
    k_gas(id_T) = 0.05;

   
    end
    % Fit the data to the polynomial and store the coefficients in the matrix:
    matrixProperties{1, id_compound}=Polyfit(T_liq, rho_liq, poly_order);
    matrixProperties{2, id_compound}=Polyfit(T_liq, C_liq, poly_order);
    matrixProperties{3, id_compound}=Polyfit(T_liq, k_liq, poly_order);
    matrixProperties{4, id_compound}=Polyfit(T_liq, mu_liq, poly_order);
    matrixProperties{5, id_compound}=Polyfit(T_gas, rho_gas, poly_order);
    matrixProperties{6, id_compound}=Polyfit(T_gas, Cp_gas, poly_order);
    matrixProperties{7, id_compound}=Polyfit(T_gas, k_gas, poly_order);

    
    % Optional: Check Polyfit accuracy (disabled by default)
    if verification_plots==1
        figure(1)
        plot(T_liq, rho_liq, 'Marker', 'o')
        hold on
        plot(T_liq, Polyval(Polyfit(T_liq, rho_liq, poly_order), T_liq))
        title(['rho liq adjustement, order: ' num2str(poly_order)])

        figure(2)
        plot(T_liq, C_liq, 'Marker', 'o')
        hold on
        plot(T_liq, Polyval(Polyfit(T_liq, C_liq, poly_order), T_liq))
        title(['C liq adjustement, order: ' num2str(poly_order)])

        figure(3)
        plot(T_liq, k_liq, 'Marker', 'o')
        hold on
        plot(T_liq, Polyval(Polyfit(T_liq, k_liq, poly_order), T_liq))
        title(['k liq adjustement, order: ' num2str(poly_order)])

        figure(4)
        plot(T_liq, mu_liq, 'Marker', 'o')
        hold on
        plot(T_liq, Polyval(Polyfit(T_liq, mu_liq, poly_order), T_liq))
        title(['mu liq adjustement, order: ' num2str(poly_order)])

        figure(5)
        plot(T_gas, rho_gas, 'Marker', 'o')
        hold on
        plot(T_gas, Polyval(Polyfit(T_gas, rho_gas, poly_order), T_gas))

        figure(6)
        plot(T_gas, Cp_gas, 'Marker', 'o')
        hold on
        plot(T_gas, Polyval(Polyfit(T_gas, Cp_gas, poly_order), T_gas))

        figure(7)
        plot(T_gas, k_gas, 'Marker', 'o')
        hold on
        plot(T_gas, Polyval(Polyfit(T_gas, k_gas, poly_order), T_gas))


    end

end

if verification_plots==1
    % Add legends for each series:
    figure(1) %rho_liq: only fuels
    legend(repelem([gota.tipo_combustible],2))

    figure(2) %C_liq: only fuels
    legend(repelem([gota.tipo_combustible],2))

    figure(3) %k_liq: only fuels
    legend(repelem([gota.tipo_combustible],2))
    
    figure(4) %mu_liq: only fuels
    legend(repelem([gota.tipo_combustible],2))

    figure(5) %rho_gas: both inerts and fuels
    legend(repelem([comp_inerts, gota.tipo_combustible],2))

    figure(6) %Cp_gas: both inerts and fuels
    legend(repelem([comp_inerts, gota.tipo_combustible],2))

    figure(7) %k_gas: both inerts and fuels
    legend(repelem([comp_inerts, gota.tipo_combustible],2))

end 

% toc

end