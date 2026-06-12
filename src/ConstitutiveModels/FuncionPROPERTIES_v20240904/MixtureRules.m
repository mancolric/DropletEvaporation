function retArray = MixtureRules(property, T_eval, Yi_eval, Yf_eval, MatrixProperties, gota, comp_inerts, Pressure)
% Function that estimates any target property for any MIXTURE (Y_eval) at
% any temperature (T_eval).

N_inerts=size(comp_inerts,2);
N_fuels=size(gota.tipo_combustible,2);

% (initial checkup: Yf+Yi=1 for all the nodes)
% for id_node=1:size(T_eval,2)
%     if abs(sum(Yi_eval(:,id_node)) + sum(Yf_eval(:,id_node))-1) > 1e-6
%         disp('ERROR: Mass fractions not equal to 1!')
%         return
%     end
% end
if any(abs(sum(Yi_eval, 1) + sum(Yf_eval, 1) - 1) > 1e-6)
    disp('ERROR: Mass fractions not equal to 1!')
    return
end


% Calculate the property for the target mixture at the target temperature:
switch property
    
    case 'rho_liq'
        % First obtain the individual rho_liq values at the target temperature
        % and multiply them by the corresponding mass fraction:
        row_rho_liq=1;
        rho_liq_individual=zeros(N_inerts+N_fuels, size(T_eval,2));
        rho_liq_individual_ponderated=zeros(N_inerts+N_fuels, size(T_eval,2));
        for id_i=N_inerts+1:N_inerts+N_fuels
            polynomial=squeeze(MatrixProperties{row_rho_liq, id_i});
            % Tb = gota.Tb(id_i-N_inerts);
            % T_eval_act = T_eval;
            % T_eval_act(T_eval > Tb) = Tb;
            rho_liq_individual(id_i, :) = Polyval(polynomial, T_eval);
            rho_liq_individual_ponderated(id_i, :) = Yf_eval(id_i-N_inerts, :)./rho_liq_individual(id_i, :);
        end
        % According to (Poling et al., 2001) the rho_liq of the mixture is simply:
        retArray = 1./sum(rho_liq_individual_ponderated);
        
        
    case 'C_liq'
        % First obtain the individual C_liq values at the target temperature
        % and multiply them by the corresponding mass fraction:
        row_C_liq=2;
        C_liq_individual=zeros(N_inerts+N_fuels, size(T_eval,2));
        C_liq_individual_ponderated=zeros(N_inerts+N_fuels, size(T_eval,2));
        for id_i=N_inerts+1:N_inerts+N_fuels
            polynomial=squeeze(MatrixProperties{row_C_liq, id_i});
            % Tb = gota.Tb(id_i-N_inerts);
            % T_eval_act = T_eval;
            % T_eval_act(T_eval > Tb) = Tb;
            C_liq_individual(id_i, :) = Polyval(polynomial, T_eval);
            C_liq_individual_ponderated(id_i, :) = C_liq_individual(id_i, :).*Yf_eval(id_i-N_inerts, :);
        end
        % The Cp_gas of the mixture is simply the sum of all columns (i.e., of all compounds):
        retArray = sum(C_liq_individual_ponderated);
        
        
    case 'k_liq'
        row_k_liq=3;
        row_rho_liq=1;
        % Li method (recommended in (Perry and Green, 2008) for mixtures above 2 compounds; implemented here due to its general applicability):
        
        % Step 1: obtain the Xf_eval, rho_liq and k_liq:
        [~, Xf_eval] = calculate_moleFractions(comp_inerts,Yi_eval, gota, Yf_eval);
        rho_liq_individual=zeros(N_inerts+N_fuels, size(T_eval,2));
        k_liq_individual=zeros(N_inerts+N_fuels, size(T_eval,2));
        for id_f=N_inerts+1:N_inerts+N_fuels
            polynomial_rholiq=squeeze(MatrixProperties{row_rho_liq, id_f});
            polynomial_kliq=squeeze(MatrixProperties{row_k_liq, id_f});
            % Tb = gota.Tb(id_f-N_inerts);
            % T_eval_act = T_eval;
            % T_eval_act(T_eval > Tb) = Tb;
            rho_liq_individual(id_f, :) = Polyval(polynomial_rholiq, T_eval);
            k_liq_individual(id_f, :) = Polyval(polynomial_kliq, T_eval);
        end
        
        % Step 2: Li averaging. Due to its complexity, it is calculated point by
        % point:
        % retArray = zeros(1, size(Xf_eval,2));
        % for id_node=1:size(Xf_eval,2)
        %     retArray(id_node) = calcula_lambda_liq_mix(gota, Xf_eval(:,id_node), rho_liq_individual(N_inerts+1:N_inerts+N_fuels,id_node), k_liq_individual(N_inerts+1:N_inerts+N_fuels,id_node));
        % end
        retArray = calcula_lambda_liq_mix(gota, Xf_eval, rho_liq_individual(N_inerts+1:N_inerts+N_fuels,:), k_liq_individual(N_inerts+1:N_inerts+N_fuels,:));
        
    case 'D_liq'
        % Wilke-Chang approximation proposed in (Sazhin et al. (2014)) as a
        % compromise between accuracy and simplicity. Applicable to any number
        % of compounds. This approximate method yields the same diffusion
        % coefficients for all the compounds in the mixture.
        
        % Step 1: obtain mu_liq through the Grunberg and Nissan rule (recommended in Poling et al. (2001)),
        % Molar fractions required:
        [~, Xf_eval] = calculate_moleFractions(comp_inerts,Yi_eval, gota, Yf_eval);
        
        row_mu_liq=4;
        mu_liq_individual=zeros(N_inerts+N_fuels, size(T_eval,2));
        mu_liq_individual_ponderated=zeros(N_inerts+N_fuels, size(T_eval,2));
        for id_i=N_inerts+1:N_inerts+N_fuels
            polynomial=squeeze(MatrixProperties{row_mu_liq, id_i});
            % Tb = gota.Tb(id_i-N_inerts);
            % T_eval_act = T_eval;
            % T_eval_act(T_eval > Tb) = Tb;
            mu_liq_individual(id_i, :) = Polyval(polynomial, T_eval);
            mu_liq_individual_ponderated(id_i, :) = log(mu_liq_individual(id_i, :)).*Xf_eval(id_i-N_inerts, :);
        end
        mu_liq_mix = exp(sum(mu_liq_individual_ponderated));
        
        % Step 2: obtain averaged values for MW and V:
        MW_avg=1000*sum(gota.mw'.*Xf_eval); % Average molar mass (kg/kmol)
        V_v = ((1.468*MW_avg.^0.297)/1.18).^3; % Molecular volume. Eq. 40 paper (Sazhin et al. (2014))
        
        % Step 3: obtain D_liq:
        retArray = 7.4e-15*T_eval.*(MW_avg).^0.5./(mu_liq_mix.*V_v.^0.6); %(m2/s)
        
        
    case 'rho_gas'
        % Two options:
        % a) Mixture rule on the fitted polynomials available in
        % MatrixProperties (better suited for real gases)
        
        % b) For ideal gas: direct calculation through ideal gas equation:
        [Xi_eval, Xf_eval] = calculate_moleFractions(comp_inerts,Yi_eval, gota, Yf_eval);
        props_kinetic_theory=zeros(N_inerts, 5);
        MW_i = zeros(N_inerts, 1);
        for id_i=1:N_inerts
            props_kinetic_theory(id_i,:) = propiedades_Tcinetica_gasMonocomp(comp_inerts(id_i));
            MW_i(id_i, 1) = props_kinetic_theory(id_i,1); %kg/mol
        end
        if N_fuels>1
            MW_avg=sum(MW_i.*Xi_eval)+sum(gota.mw'.*Xf_eval); % Average molar mass (kg/mol)
        elseif N_fuels==1
            MW_avg=sum(MW_i.*Xi_eval)+gota.mw'.*Xf_eval; % Average molar mass (kg/mol)
        end
        Ru=8.314; % J/mol/K
        retArray = Pressure*MW_avg./(Ru*T_eval);
        
        %%%%%%%%%%%%%%%%%%%% Cambiado %%%%%%%%%%%%%%%%%%%%
        
        % retArray = 0.5*ones(size(retArray));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    case 'Cp_gas'
        % First: obtain the individual Cp values at the target temperature
        % and multiply them by the corresponding mass fraction:
        row_Cp_gas=6;
        Cp_gas_individual=zeros(N_inerts+N_fuels, size(T_eval,2));
        Cp_gas_individual_ponderated=zeros(N_inerts+N_fuels, size(T_eval,2));
        for id_i=1:N_inerts
            polynomial=squeeze(MatrixProperties{row_Cp_gas, id_i});
            Cp_gas_individual(id_i, :) = Polyval(polynomial, T_eval);
            Cp_gas_individual_ponderated(id_i, :) = Cp_gas_individual(id_i, :).*Yi_eval(id_i, :);
        end
        for id_i=N_inerts+1:N_inerts+N_fuels
            polynomial=squeeze(MatrixProperties{row_Cp_gas, id_i});
            Cp_gas_individual(id_i, :) = Polyval(polynomial, T_eval);
            Cp_gas_individual_ponderated(id_i, :) = Cp_gas_individual(id_i, :).*Yf_eval(id_i-N_inerts, :);
        end
        % The Cp_gas of the mixture is simply the sum of all columns (i.e., of all compounds):
        retArray = sum(Cp_gas_individual_ponderated);
        
        
    case 'k_gas'
        % Molar fractions required:
        [Xi_eval, Xf_eval] = calculate_moleFractions(comp_inerts,Yi_eval, gota, Yf_eval);
        
        % Mixture rule applied: 'Mixture-averaged formula from Kee' (Kee et al. (2005))
        row_k_gas=7;
        k_gas_individual=zeros(N_inerts+N_fuels, size(T_eval,2));
        term_1=zeros(N_inerts+N_fuels, size(T_eval,2));
        term_2=zeros(N_inerts+N_fuels, size(T_eval,2));
        for id_i=1:N_inerts
            polynomial=squeeze(MatrixProperties{row_k_gas, id_i});
            k_gas_individual(id_i, :) = Polyval(polynomial, T_eval);
            term_1(id_i, :) =  k_gas_individual(id_i, :).*Xi_eval(id_i, :); % X_i*k_i
            term_2(id_i, :) =  Xi_eval(id_i, :)./k_gas_individual(id_i, :); % X_i/k_i
        end
        for id_i=N_inerts+1:N_inerts+N_fuels
            polynomial=squeeze(MatrixProperties{row_k_gas, id_i});
            k_gas_individual(id_i, :) = Polyval(polynomial, T_eval);
            term_1(id_i, :) =  k_gas_individual(id_i, :).*Xf_eval(id_i-N_inerts, :);
            term_2(id_i, :) =  Xf_eval(id_i-N_inerts, :)./k_gas_individual(id_i, :); % X_i/k_i
        end
        % The k_gas of the mixture is:
        retArray = 0.5*(sum(term_1)+1./(sum(term_2)));
        
    case 'ap_gas'
        % First: obtain the individual ap values at the target
        % temperature:
        row_ap_gas=8;
        ap_gas_individual  = zeros(N_inerts+N_fuels, size(T_eval,2));
        retArray           = cell(N_inerts+N_fuels,1);
        for id_i=1:N_inerts
            polynomial                  =squeeze(MatrixProperties{row_ap_gas, id_i});
            ap_gas_individual(id_i, :)  = Polyval(polynomial, T_eval);
            retArray{id_i}              = ap_gas_individual(id_i, :);
        end
        for id_i=N_inerts+1:N_inerts+N_fuels
            polynomial                  =squeeze(MatrixProperties{row_ap_gas, id_i});
            ap_gas_individual(id_i, :)  = Polyval(polynomial, T_eval);
            retArray{id_i}              = ap_gas_individual(id_i, :);
        end
        
    case 'D_gas'
        % Step 1: obtain the binary diffusion coefficients between fuel
        % (averaged kinetic theory properties) and the 4 inert compounds.
        % Method: kinetic theory of gases detailed in (Kee et al. (2005)).
        retArray           = cell(N_inerts+N_fuels,1);
        
        [Xi_eval, Xf_eval] = calculate_moleFractions(comp_inerts,Yi_eval, gota, Yf_eval);
        Xf_eval            = max(Xf_eval, 1e-12);
        if N_fuels>1
            epsilon_f_eval = Xf_eval./sum(Xf_eval);
        elseif N_fuels==1
            epsilon_f_eval = Xf_eval./Xf_eval;
        end
        
        D_fN2       = calcula_D_binario_Kee_vector2(gota, 'N2', T_eval, epsilon_f_eval);
        D_fO2       = calcula_D_binario_Kee_vector2(gota, 'O2', T_eval, epsilon_f_eval);
        D_fCO2      = calcula_D_binario_Kee_vector2(gota, 'CO2', T_eval, epsilon_f_eval);
        D_fH2O      = calcula_D_binario_Kee_vector2(gota, 'H2O', T_eval, epsilon_f_eval);
        D_fCO       = calcula_D_binario_Kee_vector2(gota, 'CO', T_eval, epsilon_f_eval);
        
        D_N2O2      = calcula_D_binario_Kee_vector2('N2', 'O2', T_eval, epsilon_f_eval);
        D_N2CO2     = calcula_D_binario_Kee_vector2('N2', 'CO2', T_eval, epsilon_f_eval);
        D_N2H2O     = calcula_D_binario_Kee_vector2('N2', 'H2O', T_eval, epsilon_f_eval);
        D_N2CO      = calcula_D_binario_Kee_vector2('N2', 'CO', T_eval, epsilon_f_eval);
        
        D_O2CO2     = calcula_D_binario_Kee_vector2('O2', 'CO2', T_eval, epsilon_f_eval);
        D_O2H2O     = calcula_D_binario_Kee_vector2('O2', 'H2O', T_eval, epsilon_f_eval);
        D_O2CO      = calcula_D_binario_Kee_vector2('O2', 'CO', T_eval, epsilon_f_eval);
        
        D_CO2H2O    = calcula_D_binario_Kee_vector2('CO2', 'H2O', T_eval, epsilon_f_eval);
        D_CO2CO     = calcula_D_binario_Kee_vector2('CO2', 'CO', T_eval, epsilon_f_eval);
        
        D_H2OCO     = calcula_D_binario_Kee_vector2('H2O', 'CO', T_eval, epsilon_f_eval);
        
        % Step 2: Calculate a global mass diffusion coefficient through the Wilke approximation
        % (Fairbanks and Wilke, 1950). This simplified aproach allows for a
        % single mass diff. coeff., which is common to all species.
        if N_fuels>1
            X_f_tot=sum(Xf_eval);
            Y_f_tot=sum(Yf_eval);
        elseif N_fuels==1
            X_f_tot=Xf_eval;
            Y_f_tot=Yf_eval;
        end
        X_N2=Xi_eval(1, :);
        X_O2=Xi_eval(2, :);
        X_CO2=Xi_eval(3, :);
        X_H2O=Xi_eval(4, :);
        X_CO=Xi_eval(5, :);
        
        %%%%%%%%% Cambiado %%%%%%%%%
        Y_N2=Yi_eval(1, :);
        Y_O2=Yi_eval(2, :);
        Y_CO2=Yi_eval(3, :);
        Y_H2O=Yi_eval(4, :);
        Y_CO=Yi_eval(5, :);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        retArray{1}     = (1 - X_N2) ./(X_O2./D_N2O2  + X_CO2./D_N2CO2 + X_H2O./D_N2H2O  + X_CO./D_N2CO    + X_f_tot./D_fN2);
        retArray{2}     = (1 - X_O2) ./(X_N2./D_N2O2  + X_CO2./D_O2CO2 + X_H2O./D_O2H2O  + X_CO./D_O2CO    + X_f_tot./D_fO2);
        retArray{3}     = (1 - X_CO2)./(X_N2./D_N2CO2 + X_O2./D_O2CO2  + X_H2O./D_CO2H2O + X_CO./D_CO2CO   + X_f_tot./D_fCO2);
        retArray{4}     = (1 - X_H2O)./(X_N2./D_N2H2O + X_O2./D_O2H2O  + X_CO2./D_CO2H2O + X_CO./D_H2OCO   + X_f_tot./D_fH2O);
        retArray{5}     = (1 - X_CO)./(X_N2./D_N2CO + X_O2./D_O2CO     + X_CO2./D_CO2CO  + X_H2O./D_H2OCO  + X_f_tot./D_fCO);
        
        retArray{6}     = (1 - X_f_tot)./(X_O2./D_fO2 + X_N2./D_fN2 + X_H2O./D_fH2O + X_CO2./D_fCO2 + X_CO./D_fCO);
        
    otherwise
        error('Property not found: %s', property)
end
end