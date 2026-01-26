function value = devuelve_propiedad_fuel(obj, propiedad, comb, estado, T )
% [VERSION MODELO NUMERICO A.GUT 31/10/2025]


% Variable que activa/desactiva warnings por evaluar propiedades fuera de
% rango:
muestra_aviso=0;


% Se cargan las propiedades básicas necesarias para el compuesto en cuestión:
% Para ello hay que buscar el índice dentro de obj.tipo_combustible:
for iD_combustible=1:size(obj.tipo_combustible,2)
    if strcmpi(obj.tipo_combustible{iD_combustible}, comb)
        indice_comb=iD_combustible;
        break
    end
end
MW=obj.mw(indice_comb)*1000;    %  OJO! Aqui en [g/mol]
Tc=obj.Tc(indice_comb);         % [K]
Tb= obj.Tb(indice_comb);        % [K]
Tbr= Tb/Tc;                     % [adimensional]
% Tr=T/Tc;% [adimensional] --> Se coloca después de fijar Tmax=Tb 
Pc= obj.Pc(indice_comb)*1e-5;   % OJO! Aqui en [bar]
R= 8.314;                       % Constante de los gases ideales, en J/(mol.K)
T_fusion = obj.Tf(indice_comb);

% IMPORTANTE: EVALUACION DE PROPIEDADES FUERA DE LIMITES:
% Dejar de considerar dependencia de propiedades con T en caso de que
% T_evaluacion > T_c! (aplicable en mezclas multicomponentes)
if strcmpi(estado, 'liquido') && T>=Tc-2
    T=Tc-2;    
%     if muestra_aviso==1
%         fprintf('### [%s] WARNING: T_eval > T_b when evaluating liquid properties! T_eval has been changed to T_b ### \n ',char(comb))
%     end
elseif strcmpi(estado, 'liquido') && T<T_fusion % El liquido se congelaría!
    % Las propiedades se extrapolan, pero la simulación dejaría de ser
    % válida. Se interrumpe aquí:
    % fprintf('### [%s] WARNING: T_eval > T_freezing when evaluating liquid properties! Code will break! ### \n ',char(comb))
    % return % Se interrumpe el modelo (eliminar/comentar "return" si se quiere proseguir simulación extrapolando propiedades del liquido para T<Tf)
end


% Calculo la Temp. reducida
% (se pone después del chequeo T>Tb)
Tr=T/Tc;% [adimensional]


% #########################################################################
% Propiedades HEPTANO: Kee-CHEMKIN (lambda_v, Dif_mix), NASA Polynomials
% (Cpv), NIST (Lv, Pvap), Perrys (resto), NASA (Hf)
% #########################################################################
if strcmpi(comb,'heptano')
    
    w= 0.350;   % Factor acéntrico
    dipole= 0;  % Debye, [D]
    Vc= 428;    % [cm3/mol]
    
    % Propiedades HEPTANO(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        if strcmpi(propiedad,'Pvap')
            % Presión de vapor del líquido, 'pvap' [Pa]
            %P_vapor, NIST
            if (290 <= T)&&(T <= 542)
                % A=4.02832; B=1268.636; C=-56.199;
                % value = 10^5*10^(A-B/(T+C)); %(Pa)
                c1=-7.897398; c1_5=2.866; c2=-2.990959; c3=1.196037; c4=-4.925711;
                Pc=2734300;
                Tc=540.13;
                tau = 1 - T/Tc;
                value = exp(log(Pc) + (Tc/T) * (c1*tau + c1_5*tau^1.5 + c2*tau^2 + c3*tau^3 + c4*tau^4)); %(Pa)
            else % Compruebo que la extrapolación es razonable...
                A=4.02832; B=1268.636; C=-56.199;
                value = 10^5*10^(A-B/(T+C)); %(Pa)
                if muestra_aviso==1
                    disp('### [Heptano] AVISO: T fuera de rango en cálculo Pvap (< 298K). Se extrapolará ###')
                end
            end
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
        elseif strcmpi(propiedad,'Lv')
            % Lv, NIST
            A=53.66; beta=0.2831;
            if (298 <= T)&&(T <= 510)
                Tr=T/Tc;
                value = 1e6*(A*exp(-beta*Tr)*(1-Tr)^beta)/MW; %J/kg
            else % Compruebo que la extrapolación es razonable...
                Tr=T/Tc;
                value = 1e6*(A*exp(-beta*Tr)*(1-Tr)^beta)/MW; %J/kg
                if muestra_aviso==1
                    disp('### [Heptano] AVISO: T fuera de rango en cálculo Lv. Se extrapolará(< 298K) ###')
                end
            end
            
            % Densidad del líqido, 'rho' [kg/m3]
        elseif strcmpi(propiedad,'rho')
            % rho_l [kg/m3]
            % Perrys Handbook
            C1=0.61259; C2=0.26211; C3=540.2;C4=0.28141;
            if (182.57 <= T)&&(T <= 550)
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
            else 
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
                % Poco probable que T<182 K... 
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            
            % Viscosidad, 'mu' [Pa.s]
        elseif strcmpi(propiedad,'mu')
            % Perrys Handbook
            C1=-9.4622; C2=877.07; C3=-0.23445; C4=1.4022E+22; C5=-10;
            if (180.15 <= T)&&(T <= 510)
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else 
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
                % Poco probable que T<182 K..
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Calor específico del liquido, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            % Perrys Handbook
            C1=61.26; C2=314410; C3=1824.6; C4=-2547.9; C5=0; a=1-T/Tc;
            if (182.57 <= T)&&(T <= 510)
                % Nota: Para el heptano, ec. 2, para el resto, ec. 1.
                value=((C1^2)/a + C2 - 2*C1*C3*a - C1*C4*a^2 - (C3^2)*(a^3)/3 - C3*C4*(a^4)/2 - (C4^2)*(a^5)/5)/MW;
            else
                % AJUSTE A POLINOMIO HECHO POR A.GUT. OCTUBRE 2025
                C1 = 18313143.336457685; C2 = -72997.88453983126; C3 = 74.14574127249222;
                value = (C1 + C2*T + C3*T^2) / MW;
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Perrys Handbook. Relacion lineal con T, por lo que no parece
            % peligroso extrapolar...
            C1=0.215; C2=-0.000303; C3=0; C4=0; C5=0;
            value=C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4;

        % Entalpía de formación, 'Hf' [J/Kg]
        % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
                value=-187.78*1000/(MW/1000);
        end
        
        
        % Propiedades HEPTANO(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        if strcmpi(propiedad,'mu')
            % Perrys Handbook
            C1=6.672E-08; C2=0.82837; C3= 85.752; C4= 0;
            if (182.57 <= T)&&(T <= 1000)
                value=C1*(T^C2)/(1+C3/T+C4/(T^2));
            elseif T>1000 % Method of Chung et al., Poling et al., The Properties of Gases
                % and Liquids, pp. 9.7.
                % [Al  cambiar de un metodo a otro hay una discontinuidad, pero la mu_v tiene una influencia prácticamente nula en el modelo (calculo Re)]
                reduced_dipole=131.3*dipole/((Vc*Tc)^(1/2));
                Fc=1-0.2756*w+0.059035*reduced_dipole^4;
                T_ast=1.2593*Tr;
                A=1.16145; B=0.14874; C=0.52487; D=0.77320; E=2.16178; F=2.43787;
                Omega_v=(A*T_ast^(-B))+C*exp(-D*T_ast)+E*exp(-F*T_ast);
                if (0.3 <= T_ast)&&(T_ast <= 100)
                    muG=40.785*Fc*((MW*T)^(1/2))/(Omega_v*Vc^(2/3)); % Micropoise
                    value=muG*10^(-7); % Pa.s
                    % disp('### AVISO: Método de estimación en curso ###')
                else
                    disp('### AVISO: Temperatura fuera de rango ###')
                end
            end

            % Entalpía de formación, 'Hf' [J/Kg]
        elseif strcmpi(propiedad,'Hf')
                value=-187.78*1000/(MW/1000);

            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            % NASA polynomial fit
            if (200 <= T)&&(T <= 1000) %Primer intervalo
                a1=1.11532484E+01; a2=-9.49415433e-3; a3=1.95571181E-04; a4=-2.49752520E-07; a5=9.84873213E-11;
                value=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/MW;
            elseif (T > 1000) %Segundo intervalo
                a1=1.85354704e+01; a2=3.91420468e-2; a3=-1.38030268e-5; a4=2.22403874E-09; a5=-1.33452580e-13;
                value=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/MW;
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Conductividad vapor fuel (W/mK): Kee-CHEMKIN. Teoria cinética
            % gases (no tienen limite de T):
            value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
        end
    end
    
% #########################################################################
% Propiedades HEXANO:  Kee-CHEMKIN (lambda_v, Dif_mix), 
%  NIST (Lv, Pvap), Perrys (resto), NASA (Hf)
% #########################################################################
elseif strcmpi(comb,'hexano')
      
    % Propiedades HEXANO(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        if strcmpi(propiedad,'Pvap')
            % Presión de vapor del líquido, 'pvap' [Pa]
            %P_vapor, NIST
            if (286 <= T)&&(T <= Tb)
                A=4.00266; B=1171.53; C=-48.784;
                value = 10^5*10^(A-B/(T+C)); %(Pa)
            else % Compruebo que la extrapolación es razonable...
                A=4.00266; B=1171.53; C=-48.784;
                value = 10^5*10^(A-B/(T+C)); %(Pa)
                if muestra_aviso==1
                    disp('### [Hexano] AVISO: T fuera de rango en cálculo Pvap. Se extrapolará ###')
                end
            end
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
        elseif strcmpi(propiedad,'Lv')
            % Lv, NIST
            A=43.85; alfa=-0.039; beta=0.397;
            if (298 <= T)&&(T <= 444)
                value = 1e6*(A*exp(-alfa*Tr)*(1-Tr)^beta)/MW; %J/kg
            else  % Compruebo que la extrapolación es razonable...
                value = 1e6*(A*exp(-alfa*Tr)*(1-Tr)^beta)/MW; %J/kg
                if muestra_aviso==1
                    disp('### [Hexano] AVISO: T fuera de rango en cálculo Lv. Se extrapolará ###')
                end
            end
            
            % Densidad del líqido, 'rho' [kg/m3]
        elseif strcmpi(propiedad,'rho')
            % rho_l [kg/m3]
            % Perrys Handbook
            C1=0.70824; C2=0.26411; C3=507.6; C4=0.27537;
            if (177.8 <= T)&&(T <= Tb)
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
            else 
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
                % Poco probable que T<177 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Viscosidad, 'mu' [Pa.s]
        elseif strcmpi(propiedad,'mu')
            % Perrys Handbook
            C1=-6.3276; C2=640; C3=-0.694; C4=5.6884E+21; C5=-10;
            if (174 <= T)&&(T <= Tb)
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else 
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
                % Poco probable que T<174 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Calor específico del liquido, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            % Perrys Handbook
            C1=172.120; C2=-183.78; C3=0.88734; C4=0; C5=0;
            if (178 <= T)&&(T <= Tb)
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
            else 
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
                % Poco probable que T<178 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Perrys Handbook. Relacion lineal con T, por lo que no parece
            % peligroso extrapolar...
            C1=0.22492; C2=-0.0003533; C3=0; C4=0; C5=0;
            value=C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4;

            % Entalpía de formación, 'Hf' [J/Kg]
            % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
            value=-166.92*1000/(MW/1000);
        end
        
        
        
        % Propiedades HEXANO(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        if strcmpi(propiedad,'mu')
            % Perrys Handbook
            C1=1.7514E-07; C2=0.70737; C3=157.14; C4= 0;
            if (178 <= T)&&(T <= 1000)
                value=C1*(T^C2)/(1+C3/T+C4/(T^2));
            else % Se comprueba que la extrapolacion es razonable...
                value=C1*(T^C2)/(1+C3/T+C4/(T^2));
                if muestra_aviso==1
                    disp('### [Hexano] AVISO: T fuera de rango en cálculo mu_v. Se extrapolará ###')
                end
            end
            
            % Entalpía de formación, 'Hf' [J/Kg]
            % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
            value=-166.92*1000/(MW/1000);
            
          % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp') %PERRYS (no está en el NASA Polynomials)
            C1=1.044E+05; C2=3.523E+05; C3=1.6946E+03; C4=2.369E+05; C5=761.6;
            if (200 <= T)&&(T <= 1500)
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
            else % Se comprueba que la extrapolacion es razonable...
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
                if muestra_aviso==1
                    disp('### [Hexano] AVISO: T fuera de rango en cálculo Cp_v. Se extrapolará ###')
                end
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Conductividad vapor fuel (W/mK): Kee-CHEMKIN. Teoria cinética
            % gases (no tienen limite de T):
            value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
        end
    end
    
    % #########################################################################
    % ______Propiedades OCTANO: Kee-CHEMKIN (lambda_v, Dif_mix), Perrys (resto)
    %                           NASA (Hf)
    % #########################################################################
elseif strcmpi(comb,'octano') 
    
    w= 0.399; % Factor acéntrico
    dipole= 0;  % Debye, [D]
    Vc= 492;    % [cm3/mol]
    
    % Propiedades OCTANO(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        % Presión de vapor del líquido, 'pvap' [Pa]
        if strcmpi(propiedad,'pvap') % Perrys Handbook (NIST lo da a tramos...)
            C1=96.084; C2=-7900.2; C3=-11.003; C4=7.1802E-06; C5=2;
            if (216.38<=T)&&(T<=568.7)
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else 
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
                % Poco probable que T<216 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
        elseif strcmpi(propiedad,'Lv') % Perrys Handbook (NIST lo da a tramos...)
            C1=5.518E+07; C2=0.38467; C3=0; C4=0; C5=0;
            if (216.38 <= T)&&(T <= 568.7)
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
            else
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
                % Poco probable que T<216 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Densidad del líqido, 'rho' [kg/m3]
        elseif strcmpi(propiedad,'rho')  % Perrys Handbook (
            C1=0.5266; C2=0.25693; C3=568.7; C4=0.28571;
            if (216.38 <= T)&&(T <= 568.7)
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
            else
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
                % Poco probable que T<216 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Viscosidad, 'mu' [Pa.s]
        elseif strcmpi(propiedad,'mu')
            C1=-7.556; C2=881.09; C3=-0.52502; C4=4.6342E+22; C5=-10;
            if (211.15 <= T)&&(T <= 454.96)
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
                % Poco probable que T<216 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Calor específico, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=224830; C2=-186.63; C3=0.95891; C4=0; C5=0;
            if (216.38 <= T)&&(T <= 460)
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
            else
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
                % Poco probable que T<216 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Perrys Handbook. Relacion lineal con T, por lo que no parece
            % peligroso extrapolar...
            C1=0.2156; C2=-0.00029483; C3=0; C4=0; C5=0;
            value=C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4;

            % Entalpía de formación, 'Hf' [J/Kg]
            % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
            value=-208.75*1000/(MW/1000);
        end
        
        
        % Propiedades OCTANO(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        if strcmpi(propiedad,'mu')
            C1=3.1191E-08; C2=0.92925; C3= 55.092; C4= 0;
            if (216.38 <= T)&&(T <= 1000)
                value=C1*(T^C2)/(1+C3/T+C4/(T^2));
            elseif T>1000 % Method of Chung et al., Poling et al., The Properties of Gases
                % and Liquids, pp. 9.7
                reduced_dipole=131.3*dipole/((Vc*Tc)^(1/2));
                Fc=1-0.2756*w+0.059035*reduced_dipole^4;
                T_ast=1.2593*Tr;
                A=1.16145; B=0.14874; C=0.52487; D=0.77320; E=2.16178; F=2.43787;
                Omega_v=(A*T_ast^(-B))+C*exp(-D*T_ast)+E*exp(-F*T_ast);
                if (0.3 <= T_ast)&&(T_ast <= 100)
                    muG=40.785*Fc*((MW*T)^(1/2))/(Omega_v*Vc^(2/3)); % Micropoise
                    value=muG*10^(-7); % Pa.s
                    % disp('### AVISO: Método de estimación en curso ###')
                else
                    disp('### AVISO: Temperatura fuera de rango ###')
                end
            end
        elseif strcmpi(propiedad,'Hf')
            value=-208.75*1000/(MW/1000); 
            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp') % Perrys Handbook 
            C1=1.3554E+05; C2=4.431E+05; C3=1.6356E+03; C4=3.054E+05; C5=746.4;
            if (200 <= T)&&(T <= 1500)
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
            else % Compruebo extrapolacion razonable...
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
                if muestra_aviso==1
                    disp('### [Octano] AVISO: T fuera de rango en cálculo Cp_v. Se extrapolará ###')
                end
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
          % Conductividad vapor fuel (W/mK): Kee-CHEMKIN. Teoria cinética
          % gases (no tienen limite de T):
          value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
        end
    end
    
    
    
    % #########################################################################
    % _____Propiedades DODECANO: Kee-CHEMKIN (lambda_v, Dif_mix), Perrys (resto)
    % #########################################################################
elseif strcmpi(comb,'dodecano')  
    
    w= 0.576; % Factor acéntrico
    dipole= 0;  % Debye, [D]
    Vc= 754;    % [cm3/mol]
    
    % Propiedades DODECANO(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        % Presión de vapor del líquido, 'pvap' [Pa]
        if strcmpi(propiedad,'pvap') % Perrys
            C1=137.47; C2=-11976; C3=-16.698; C4=8.0906E-06; C5=2;
            if (263.57<=T)&&(T<=658)
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
                % Poco probable que T<263 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
        elseif strcmpi(propiedad,'Lv') % Perrys
            C1=7.7337E+07; C2=0.40681; C3=0; C4=0; C5=0;
            if (263.57 <= T)&&(T <= 658)
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
            else
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
                % Poco probable que T<263 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Densidad del líqido, 'rho' [kg/m3]
        elseif strcmpi(propiedad,'rho') % Perrys
            C1=0.33267; C2=0.24664; C3=658; C4=0.28571;
            if (263.57 <= T)&&(T <= 658)
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
            else
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
                % Poco probable que T<263 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Viscosidad, 'mu' [Pa.s]
        elseif strcmpi(propiedad,'mu') % Perrys
            C1=-7.8244; C2=1191.9; C3=-0.49963; C4=3.9572E+23; C5=-10;
            if (262.15 <= T)&&(T <= 526.40) 
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
                % Poco probable que T<263 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Calor específico, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=508210; C2=-1368.7; C3=3.1015; C4=0; C5=0;
            if (263.57 <= T)&&(T <= 330) % Perrys. Extrañamente tiene un limite superior de T muy bajo! En NIST no dan expresion...
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
            else % Se comprueba que la extrapolacion es razonable...
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
                if muestra_aviso==1
                    disp('### [Dodecano] AVISO: T fuera de rango en cálculo Cp_l. Se extrapolará ###')
                end
                 
             % [2020-02-13] ALTERNATIVA A EXTRAPOLAR: Corresponding states method
             % Da resultados muy parecidos a la extrapolación (lo cual confirma la validez de extrapolar),
             % pero, como es lógico, cambiar de metodo en 330 K da una
             % discontinuidad en Cp_l en 330 K. Para evitar esta
             % discontinuidad se prefiere extrapolar.
             
%             elseif (T > 330)&&(T <= 0.99*Tc) % Como mal menor, se 
%                 % Corresponding-states method, Bondi, refitted, Poling et al, 
%                 %The Properties of Gases and Liquids, pp. 6.22:
%                 
%                 % Calculamos en primer lugar el calor específico a presión
%                 % constante del vapor (lo necesitamos en J/mol.K)
%                 C1=2.1295E+05; C2=6.633E+05; C3=1.7155E+03; C4=4.5161E+05; C5=777.5;
%                 if (200 <= T)&&(T <= 1500)
%                     CpG=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/1000;
%                 else
%                     disp('### AVISO: Temperatura fuera de rango ###')
%                 end
%                 value=(CpG+R*(1.586+0.49/(1-Tr)+w*(4.2775+(6.3/Tr)*(1-Tr)^(1/3)+0.4355/(1-Tr))))*(1000/MW);
%                 % disp('### AVISO: Método de estimación en curso ###')
%             else
%                 % Poco probable que T<263 K...
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Perrys Handbook. Relacion lineal con T, por lo que no parece
            % peligroso extrapolar...
            C1=0.2047; C2=-0.0002326; C3=0; C4=0; C5=0;
            value=C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4;

            % Entalpía de formación, 'Hf' [J/Kg]
            % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
            value=-290.9*1000/(MW/1000);
        end
        
        
        % Propiedades DODECANO(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        if strcmpi(propiedad,'mu')
            C1=6.344E-08; C2=0.8287; C3= 219.5; C4= 0;
            if (263.57 <= T)&&(T <= 1000)
                value=C1*(T^C2)/(1+C3/T+C4/(T^2));
            elseif T>1000 % Method of Chung et al., Poling et al., The Properties of Gases
                % and Liquids, pp. 9.7
                reduced_dipole=131.3*dipole/((Vc*Tc)^(1/2));
                Fc=1-0.2756*w+0.059035*reduced_dipole^4;
                T_ast=1.2593*Tr;
                A=1.16145; B=0.14874; C=0.52487; D=0.77320; E=2.16178; F=2.43787;
                Omega_v=(A*T_ast^(-B))+C*exp(-D*T_ast)+E*exp(-F*T_ast);
                if (0.3 <= T_ast)&&(T_ast <= 100)
                    muG=40.785*Fc*((MW*T)^(1/2))/(Omega_v*Vc^(2/3)); % Micropoise
                    value=muG*10^(-7); % Pa.s
                    % disp('### AVISO: Método de estimación en curso ###')
                else
                    disp('### AVISO: Temperatura fuera de rango ###')
                end
            end

        elseif strcmpi(propiedad,'Hf')
            value=-290.9*1000/(MW/1000);

            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=2.1295E+05; C2=6.633E+05; C3=1.7155E+03; C4=4.5161E+05; C5=777.5;
            if (200 <= T)&&(T <= 1500)
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
            else % Compruebo extrapolacion razonable...
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
                if muestra_aviso==1
                    disp('### [Dodecano] AVISO: T fuera de rango en cálculo Cp_v. Se extrapolará ###')
                end       
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
           % Conductividad vapor fuel (W/mK): Kee-CHEMKIN. Teoria cinética
          % gases (no tienen limite de T):
           value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
        end
    end
    
    
    
    % #########################################################################
    % ___Propiedades TETRADECANO: Kee-CHEMKIN (lambda_v, Dif_mix), Perrys (resto)
    % #########################################################################
elseif strcmpi(comb,'tetradecano') 
    
    w= 0.644; % Factor acéntrico
    dipole= 0;  % Debye, [D]
    Vc= 894;    % [cm3/mol]
    
    % Propiedades TETRADECANO(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        % Presión de vapor del líquido, 'pvap' [Pa]
        if strcmpi(propiedad,'pvap')
            C1=140.47; C2=-13231; C3=-16.859; C4=6.5877E-06; C5=2;
            if (279.01<=T)&&(T<=693)
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else % Compruebo que la extrapolación es razonable...
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
                if muestra_aviso==1
                    disp('### [Tetradecano] AVISO: T fuera de rango en cálculo Pvap (< 279 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
        elseif strcmpi(propiedad,'Lv')
            C1=9.0539E+07; C2=0.44467; C3=0; C4=0; C5=0;
            if (279.01 <= T)&&(T <= 693)
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
            else % Compruebo que la extrapolación es razonable...
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
                if muestra_aviso==1
                    disp('### [Tetradecano] AVISO: T fuera de rango en cálculo Lv (< 279 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Densidad del líqido, 'rho' [kg/m3]
        elseif strcmpi(propiedad,'rho')
            C1=0.27248; C2=0.24007; C3=693; C4=0.28571;
            if (279.01 <= T)&&(T <= 693)
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
            else % Compruebo que la extrapolación es razonable...
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
                if muestra_aviso==1
                    disp('### [Tetradecano] AVISO: T fuera de rango en cálculo rho_l (< 279 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Viscosidad, 'mu' [Pa.s]
        elseif strcmpi(propiedad,'mu')
            C1=-14.493; C2=1710.8; C3=-0.4417; C4=3.0895E+28; C5=-12;
            if (277.65 <= T)&&(T <= 554.40)
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else % Compruebo que la extrapolación es razonable...
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
                if muestra_aviso==1
                    disp('### [Tetradecano] AVISO: T fuera de rango en cálculo mu_l (< 277 K = T_fus). Se extrapolará ###')
                end
                
            end
            
            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=353140; C2=29.13; C3=0.86116; C4=0; C5=0;
            if (279.01 <= T)&&(T <= 526.73)
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
            else % Compruebo que la extrapolación es razonable...
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
                if muestra_aviso==1
                    disp('### [Tetradecano] AVISO: T fuera de rango en cálculo C_l (< 279 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Perrys Handbook. Relacion lineal con T, por lo que no parece
            % peligroso extrapolar...
            C1=0.20293; C2=-0.00021798; C3=0; C4=0; C5=0;
            value=C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4;

            % Entalpía de formación, 'Hf' [J/Kg]
            % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
            value=-332.1*1000/(MW/1000);
        end
        
        
        % Propiedades TETRADECANO(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        if strcmpi(propiedad,'mu')
            C1=5.1567E-09; C2=1.1561; C3= 0; C4= 0;
            if (279.01 <= T)&&(T <= 1000)
                value=C1*(T^C2)/(1+C3/T+C4/(T^2));
            elseif T>1000 % Method of Chung et al., Poling et al., The Properties of Gases
                % and Liquids, pp. 9.7
                reduced_dipole=131.3*dipole/((Vc*Tc)^(1/2));
                Fc=1-0.2756*w+0.059035*reduced_dipole^4;
                T_ast=1.2593*Tr;
                A=1.16145; B=0.14874; C=0.52487; D=0.77320; E=2.16178; F=2.43787;
                Omega_v=(A*T_ast^(-B))+C*exp(-D*T_ast)+E*exp(-F*T_ast);
                if (0.3 <= T_ast)&&(T_ast <= 100)
                    muG=40.785*Fc*((MW*T)^(1/2))/(Omega_v*Vc^(2/3)); % Micropoise
                    value=muG*10^(-7); % Pa.s
                    % disp('### AVISO: Método de estimación en curso ###')
                else
                    disp('### AVISO: Temperatura fuera de rango ###')
                end
            end
        
        elseif strcmpi(propiedad,'Hf')
            value=-332.1*1000/(MW/1000);

            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=2.3082E+05; C2=7.8678E+05; C3=1.6823E+03; C4=5.4486E+05; C5=743.1;
            if (200 <= T)&&(T <= 1500)
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
            else % Compruebo extrapolacion razonable...
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
                if muestra_aviso==1
                    disp('### [Tetradecano] AVISO: T fuera de rango en cálculo Cp_v. Se extrapolará ###')
                end       
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
           % Conductividad vapor fuel (W/mK): Kee-CHEMKIN. Teoria cinética
          % gases (no tienen limite de T):
           value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
        end
    end
    
    
    
    % #########################################################################
    % ____Propiedades HEXADECANO: Kee-CHEMKIN (lambda_v, Dif_mix), Perrys (resto)
    % #########################################################################
elseif strcmpi(comb,'hexadecano')  
    
    w= 0.718;   % Factor acéntrico
    dipole= 0;  % Debye, [D]
    Vc= 1034;   % [cm3/mol]
    
    % Propiedades HEXADECANO(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        % Presión de vapor del líquido, 'pvap' [Pa]
        if strcmpi(propiedad,'pvap')
            C1=156.06; C2=-15015; C3=-18.941; C4=6.8172E-06; C5=2;
            if (291.31<=T)&&(T<=723)
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else % Compruebo que la extrapolación es razonable...
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
                if muestra_aviso==1
                    disp('### [Hexadecano] AVISO: T fuera de rango en cálculo Pvap (< 291 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
        elseif strcmpi(propiedad,'Lv')
            C1=10.156E+07; C2=0.45726; C3=0; C4=0; C5=0;
            if (291.31 <= T)&&(T <= 723)
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
            else % Compruebo que la extrapolación es razonable...
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
                if muestra_aviso==1
                    disp('### [Hexadecano] AVISO: T fuera de rango en cálculo Lv (< 291 K  = T_fus). Se extrapolará ###')
                end
            end
            
            % Densidad del líqido, 'rho' [kg/m3]
        elseif strcmpi(propiedad,'rho')
            C1=0.23289; C2=0.23659; C3=723; C4=0.28571;
            if (291.31 <= T)&&(T <= 723)
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
            else % Compruebo que la extrapolación es razonable...
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
                if muestra_aviso==1
                    disp('### [Hexadecano] AVISO: T fuera de rango en cálculo rho_l (< 291 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Viscosidad, 'mu' [Pa.s]
        elseif strcmpi(propiedad,'mu')
            C1=-20.182; C2=2203.5; C3=1.2289; C4=0; C5=0;
            if (291.31 <= T)&&(T <= 564.15)
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else % Compruebo que la extrapolación es razonable...
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
                if muestra_aviso==1
                    disp('### [Hexadecano] AVISO: T fuera de rango en cálculo mu_l (< 291 K = T_fus). Se extrapolará ###')
                end
                
            end
            
            % Calor específico, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=370350; C2=231.47; C3=0.68632; C4=0; C5=0;
            if (291.31 <= T)&&(T <= 560.01)
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
            else % Compruebo que la extrapolación es razonable...
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
                if muestra_aviso==1
                    disp('### [Hexadecano] AVISO: T fuera de rango en cálculo C_l (< 291 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Perrys Handbook. Relacion lineal con T, por lo que no parece
            % peligroso extrapolar...
            C1=0.20749; C2=-0.00021917; C3=0; C4=0; C5=0;
            value=C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4;

            % Entalpía de formación, 'Hf' [J/Kg]
            % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
            value=-374.9*1000/(MW/1000);

        end
        
        
        % Propiedades HEXADECANO(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        if strcmpi(propiedad,'mu')
            C1=1.2463E-07; C2=0.7322; C3=395; C4=6000;
            if (291.31 <= T)&&(T <= 1000)
                value=C1*(T^C2)/(1+C3/T+C4/(T^2));
            elseif T>1000 % Method of Chung et al., Poling et al., The Properties of Gases
                % and Liquids, pp. 9.7
                reduced_dipole=131.3*dipole/((Vc*Tc)^(1/2));
                Fc=1-0.2756*w+0.059035*reduced_dipole^4;
                T_ast=1.2593*Tr;
                A=1.16145; B=0.14874; C=0.52487; D=0.77320; E=2.16178; F=2.43787;
                Omega_v=(A*T_ast^(-B))+C*exp(-D*T_ast)+E*exp(-F*T_ast);
                if (0.3 <= T_ast)&&(T_ast <= 100)
                    muG=40.785*Fc*((MW*T)^(1/2))/(Omega_v*Vc^(2/3)); % Micropoise
                    value=muG*10^(-7); % Pa.s
                    % disp('### AVISO: Método de estimación en curso ###')
                else
                    disp('### AVISO: Temperatura fuera de rango ###')
                end
            end

        elseif strcmpi(propiedad,'Hf')
            value=-374.9*1000/(MW/1000);
            
            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=2.6283E+05; C2=8.9733E+05; C3=1.6912E+03; C4=6.264E+05; C5=744.41;
            if (200 <= T)&&(T <= 1500)
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
            else % Compruebo extrapolacion razonable...
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
                if muestra_aviso==1
                    disp('### [Hexadecano] AVISO: T fuera de rango en cálculo Cp_v. Se extrapolará ###')
                end  
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Conductividad vapor fuel (W/mK): Kee-CHEMKIN
           value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
            
        end
    end
    
    
    
    % #########################################################################
    % ____Propiedades EICOSANO: Kee-CHEMKIN (lambda_v, Dif_mix), Perrys (resto)
    % #########################################################################
elseif strcmpi(comb,'eicosano') 
    
    w= 0.865;        % Factor acéntrico
    dipole= 0;       % Debye, [D]
    %     Vref= 361.18;    % [cm3/mol]
    %     Tref= 298.15;    % [K]
    Vc= 1323.656;    % [cm3/mol]
    % NOTA: No se dispone del dato de Vc en Poling et al., The Properties
    % of Gases and Liquids. Se calcula a partir de estas dos ecuaciones:
    %   1) Zc= 0.291 - 0.080*w
    %   2) Zc= Pc*Vc/(R*Tc), siendo R= 83.145 bar.cm3/(mol.K)
    
    % Propiedades EICOSANO(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        % Presión de vapor del líquido, 'pvap' [Pa]
        if strcmpi(propiedad,'pvap')
            C1=203.66; C2=-19441; C3=-25.525; C4=8.8382E-06; C5=2;
            if (309.58<=T)&&(T<=768)
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else % Compruebo que la extrapolación es razonable...
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
                if muestra_aviso==1
                    disp('### [Eicosano] AVISO: T fuera de rango en cálculo Pvap (< 309 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
        elseif strcmpi(propiedad,'Lv')
            C1=12.86E+07; C2=0.50351; C3=0.32986; C4=-0.42184; C5=0;
            if (309.58 <= T)&&(T <= 768)
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
            else % Compruebo que la extrapolación es razonable...
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
                if muestra_aviso==1
                    disp('### [Eicosano] AVISO: T fuera de rango en cálculo Lv (< 309 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Densidad del líqido, 'rho' [kg/m3]
        elseif strcmpi(propiedad,'rho')
            C1=0.18166; C2=0.23351; C3=768; C4=0.28571;
            if (309.58 <= T)&&(T <= 768)
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
            else % Compruebo que la extrapolación es razonable...
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
                if muestra_aviso==1
                    disp('### [Eicosano] AVISO: T fuera de rango en cálculo rho_l (< 309 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Viscosidad, 'mu' [Pa.s]
        elseif strcmpi(propiedad,'mu')
            C1=-18.315; C2=2283.5; C3=0.95485; C4=0; C5=0;
            if (309.58 <= T)&&(T <= 616.93)
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else % Compruebo que la extrapolación es razonable...
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
                if muestra_aviso==1
                    disp('### [Eicosano] AVISO: T fuera de rango en cálculo mu_l (< 309 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=352720; C2=807.32; C3=0.2122; C4=0; C5=0;
            if (309.58 <= T)&&(T <= 616.93)
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
                
            else % Compruebo que la extrapolación es razonable...
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
                if muestra_aviso==1
                    disp('### [Eicosano] AVISO: T fuera de rango en cálculo C_l (< 309 K = T_fus). Se extrapolará ###')
                end

            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Perrys Handbook. Relacion lineal con T, por lo que no parece
            % peligroso extrapolar...
            C1=0.2178; C2=-0.0002233; C3=0; C4=0; C5=0;
            value=C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4;

            % Entalpía de formación, 'Hf' [J/Kg]
            % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
            value=-455.8*1000/(MW/1000);

        end
        
        % Propiedades EICOSANO(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        if strcmpi(propiedad,'mu')
            if (309.58 <= T)&&(T <= 1000)
                C1=2.9236E-07; C2=0.62458; C3=702.84; C4=0;
                value=C1*(T^C2)/(1+C3/T+C4/(T^2));
            elseif T>1000 % Method of Chung et al., Poling et al., The Properties of Gases
                % and Liquids, pp. 9.7
                reduced_dipole=131.3*dipole/((Vc*Tc)^(1/2));
                Fc=1-0.2756*w+0.059035*reduced_dipole^4;
                T_ast=1.2593*Tr;
                A=1.16145; B=0.14874; C=0.52487; D=0.77320; E=2.16178; F=2.43787;
                Omega_v=(A*T_ast^(-B))+C*exp(-D*T_ast)+E*exp(-F*T_ast);
                if (0.3 <= T_ast)&&(T_ast <= 100)
                    muG=40.785*Fc*((MW*T)^(1/2))/(Omega_v*Vc^(2/3)); % Micropoise
                    value=muG*10^(-7); % Pa.s
                    % disp('### AVISO: Método de estimación en curso ###')
                else
                    disp('### AVISO: Temperatura fuera de rango ###')
                end
            end
         
        elseif strcmpi(propiedad,'Hf')
            value=-455.8*1000/(MW/1000);

            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=3.2481E+05; C2=11.09E+05; C3=1.636E+03; C4=7.45E+05; C5=726.27;
            if (200 <= T)&&(T <= 1500)
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
            else % Compruebo extrapolacion razonable...
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
                if muestra_aviso==1
                    disp('### [Eicosano] AVISO: T fuera de rango en cálculo Cp_v. Se extrapolará ###')
                end  
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Conductividad vapor fuel (W/mK): Kee-CHEMKIN
           value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
        end
    end
    
    
    
    % #########################################################################
    % ______Propiedades METANOL: Perry's Chem. Eng. Handbook, 8th Edition._____
    % #########################################################################
elseif strcmpi(comb,'metanol')  % POR COMPLETAR
    
    %     w= 0.565; % Factor acéntrico
    %     dipole= 1.7; % Debye, [D]
    %     Vc=118;      % [cm3/mol]
    
    % Propiedades METANOL(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        % Presión de vapor del líquido, 'pvap' [Pa]
        if strcmpi(propiedad,'pvap')
            C1=82.718; C2=-6904.5; C3=-8.8622; C4=7.4664E-06; C5=2;
            if (175.47<=T)&&(T<=512.5)
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
                % Poco probable que T<175 K... 
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
        elseif strcmpi(propiedad,'Lv')
            C1=5.0451E+07; C2=0.33594; C3=0; C4=0; C5=0;
            if (175.47 <= T)&&(T <= 512.5)
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
            else
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
                % Poco probable que T<175 K... 
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Densidad del líqido, 'rho' [kg/m3]
        elseif strcmpi(propiedad,'rho')
            C1=2.3267; C2=0.27073; C3=512.5; C4=0.24713;
            if (175.47 <= T)&&(T <= 512.5)
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
            else
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
                % Poco probable que T<175 K... 
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Viscosidad, 'mu' [Pa.s]
        elseif strcmpi(propiedad,'mu')
            C1=-25.317; C2=1789.2; C3=2.069; C4=0; C5=0;
            if (175.47 <= T)&&(T <= 337.85)
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
                % Poco probable que T<175 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Calor específico liq., 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=105800; C2=-362.23; C3=0.9379; C4=0; C5=0;
            if (175.47 <= T)&&(T <= 400)
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
            else
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
                % Poco probable que T<175 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Perrys Handbook. Relacion lineal con T, por lo que no parece
            % peligroso extrapolar...
            C1=0.2837; C2=-0.000281; C3=0; C4=0; C5=0;
            value=C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4;
        
        % Entalpía de formación, 'Hf' [J/Kg]
        % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
                value=-200.94*1000/(MW/1000);
        end
        
        
        % Propiedades METANOL(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        if strcmpi(propiedad,'mu')
            C1=3.0663E-07; C2=0.69655; C3=205; C4=0;
            if (240 <= T)&&(T <= 1000) % Perrys
                value=C1*(T^C2)/(1+C3/T+C4/(T^2));
            elseif T>1000 % NASA Technical Memorandum 4513, 1993 (Expresión válida entre 1000 y 5000 K)
                A=0.60590993; B=-0.19198488E+03; C=-0.25158890E+05; D=0.17778509E+01;
                value=(10^-7)*exp(A*log(T)+B/T+C/(T^2)+D);
            end
            

        elseif strcmpi(propiedad,'Hf')
                value=-200.94*1000/(MW/1000);

            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=0.39252E+05; C2=0.879E+05; C3=1.9165E+03; C4=0.53654E+05; C5=896.7;
            if (200 <= T)&&(T <= 1500)
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
            else % Compruebo extrapolacion razonable...
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
                if muestra_aviso==1
                    disp('### [Metanol] AVISO: T fuera de rango en cálculo Cp_v. Se extrapolará ###')
                end  
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Conductividad vapor fuel (W/mK): Kee-CHEMKIN
           value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
        end
    end
    
    
    
    % #########################################################################
    % ______Propiedades ETANOL: Perry's Chem. Eng. Handbook, 8th Edition.______
    % #########################################################################
elseif strcmpi(comb,'etanol')  
    
    %     w= 0.649;    % Factor acéntrico
    %     dipole= 1.7; % Debye, [D]
    %     Vc=167;      % [cm3/mol]
    
    % Propiedades ETANOL(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        % Presión de vapor del líquido, 'pvap' [Pa]
        if strcmpi(propiedad,'pvap')
            %P_vapor, NIST
            A=5.37229; B=1670.409; C=-40.191;
            if (273 <= T)&&(T <= 352)
                value = 10^5*10^(A-B/(T+C)); %(Pa)
            else % Compruebo extrapolacion razonable...
                value = 10^5*10^(A-B/(T+C)); %(Pa)
                if muestra_aviso==1
                    disp('### [Etanol] AVISO: T fuera de rango en cálculo Pvap. Se extrapolará ###')
                end  
            end
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
            % NIST:
        elseif strcmpi(propiedad,'Lv')
            A=50.43; alfa=-0.4475; beta=0.4989;
            if (298 <= T)&&(T <= 469)
                value = 1e6*(A*exp(-alfa*Tr)*(1-Tr)^beta)/MW; %J/kg
            else % Compruebo extrapolacion razonable...
                value = 1e6*(A*exp(-alfa*Tr)*(1-Tr)^beta)/MW; %J/kg
                if muestra_aviso==1
                    disp('### [Etanol] AVISO: T fuera de rango en cálculo Lv. Se extrapolará ###')
                end  
            end
            
            % Densidad del líqido, 'rho' [kg/m3]
            % Perrys Handbook
        elseif strcmpi(propiedad,'rho')
            C1=1.6288; C2=0.27469; C3=514; C4=0.23178;
            if (159.05 <= T)&&(T <= 514)
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
            else
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
                % Poco probable que T<159 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Viscosidad, 'mu' [Pa.s]
        elseif strcmpi(propiedad,'mu')
            C1=7.875; C2=781.98; C3=-3.0418; C4=0; C5=0;
            if (200 <= T)&&(T <= 440)
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
                % Poco probable que T<200 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Calor específico liquido, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=102640; C2=-139.63; C3=-0.030341; C4=0.0020386; C5=0;
            if (159.05 <= T)&&(T <= 390)
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
            else
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
                % Poco probable que T<159 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Perrys Handbook. Relacion lineal con T, por lo que no parece
            % peligroso extrapolar...
            C1=0.2468; C2=-0.000264; C3=0; C4=0; C5=0;
            value=C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4;

        % Entalpía de formación, 'Hf' [J/Kg]
        % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
                value=-234.95*1000/(MW/1000);
        end
        
        
        % Propiedades ETANOL(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        if strcmpi(propiedad,'mu')
            C1=1.0613E-07; C2=0.8066; C3=52.7; C4=0;
            if (200 <= T)&&(T <= 1000)
                value=C1*(T^C2)/(1+C3/T+C4/(T^2));
            elseif T>1000 % NASA Technical Memorandum 4513, 1993 (Expresión válida entre 1000 y 5000 K)
                A=0.62692259; B=-0.65829493E+02; C=-0.47946575E+05; D=0.13535012E+01;
                value=(10^-7)*exp(A*log(T)+B/T+C/(T^2)+D);
            end

        elseif strcmpi(propiedad,'Hf')
                value=-234.95*1000/(MW/1000);

            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            % NASA polynomial fit
            if (200 <= T)&&(T <= 1000) %Primer intervalo
                a1=4.85868178E+00; a2=-3.74006740E-3; a3=6.95550267E-05; a4=-8.86541147E-08; a5=3.51684430E-11;
                value=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/MW;
            elseif (T > 1000) %Segundo intervalo;
                a1=6.56289770E+00; a2=1.52034264E-02; a3=-5.38922247e-6; a4=8.62150224E-10; a5=-5.12824683e-14;
                value=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/MW;
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Conductividad vapor fuel (W/mK): Kee-CHEMKIN
            value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
        end
    end
    
    
    
    % #########################################################################
    % _____Propiedades 1-BUTANOL: Perry's Chem. Eng. Handbook, 8th Edition.____
    % #########################################################################
elseif strcmpi(comb,'butanol')  
    
    w= 0.590; % Factor acéntrico
    dipole= 1.8;  % Debye, [D]
    Vc= 275;    % [cm3/mol]
    
    % Propiedades 1-BUTANOL(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        % Presión de vapor del líquido, 'pvap' [Pa]
        if strcmpi(propiedad,'pvap')
            %P_vapor, NIST
            A=4.54607; B=1351.555; C=-93.34;
            if (295.8 <= T)&&(T <= 391.0) 
                value = 10^5*10^(A-B/(T+C)); %(Pa)
            else % Compruebo extrapolacion razonable...
                value = 10^5*10^(A-B/(T+C)); %(Pa)
                if muestra_aviso==1
                    disp('### [Butanol] AVISO: T fuera de rango en cálculo Pvap. Se extrapolará ###')
                end  
            end
           
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
        elseif strcmpi(propiedad,'Lv')
             A=62.53; alfa=-0.6584; beta=0.696;
            if (298 <= T)&&(T <= 410)
                value = 1e6*(A*exp(-alfa*Tr)*(1-Tr)^beta)/MW; %J/kg
            else
                % Compruebo extrapolacion razonable...
                value = 1e6*(A*exp(-alfa*Tr)*(1-Tr)^beta)/MW; %J/kg
                if muestra_aviso==1
                    disp('### [Butanol] AVISO: T fuera de rango en cálculo Lv. Se extrapolará ###')
                end   
            end
        
            
            % Densidad del líqido, 'rho' [kg/m3]
            % Perrys Handbook
        elseif strcmpi(propiedad,'rho')
            C1=0.98279; C2=0.26830; C3=563.1; C4=0.25488;
            if (183.85 <= T)&&(T <= 563.1)
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
            else
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
                % Poco probable que T<183 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Viscosidad liquido, 'mu' [Pa.s]
        elseif strcmpi(propiedad,'mu')
            C1=0.87669; C2=1602.9; C3=-2.1475; C4=3.3866E+22; C5=-9.9231;
            if (190 <= T)&&(T <= 391.9)
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else 
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
                % Poco probable que T<190 K..
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Calor específico liquido, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=191200; C2=-730.4; C3=2.2998; C4=0; C5=0;
            if (183.85 <= T)&&(T <= 391.9)
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
            else
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
                % Poco probable que T<183 K..
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Conductividad térmica liquido, 'lambda' [W/m.K] 
        elseif strcmpi(propiedad,'lambda')
            % Perrys Handbook. Relacion lineal con T, por lo que no parece
            % peligroso extrapolar...
            C1=0.2136; C2=-0.0002034; C3=0; C4=0; C5=0;
            value=C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4;

            % Entalpía de formación, 'Hf' [J/Kg]
            % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
                value=-277.0*1000/(MW/1000);

        end
        
        
        % Propiedades 1-BUTANOL(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        if strcmpi(propiedad,'mu')
            C1=1.4031E-06; C2=0.4611; C3=537; C4=0;
            if (183.85 <= T)&&(T <= 1000)
                value=C1*(T^C2)/(1+C3/T+C4/(T^2));
            elseif T>1000 % Method of Chung et al., Poling et al., The Properties of Gases
                % and Liquids, pp. 9.7
                reduced_dipole=131.3*dipole/((Vc*Tc)^(1/2));
                kappa=0.132;
                Fc=1-0.2756*w+0.059035*reduced_dipole^4+kappa; % Término kappa: únicamente alcoholes y ácidos.
                T_ast=1.2593*Tr;
                A=1.16145; B=0.14874; C=0.52487; D=0.77320; E=2.16178; F=2.43787;
                Omega_v=(A*T_ast^(-B))+C*exp(-D*T_ast)+E*exp(-F*T_ast);
                if (0.3 <= T_ast)&&(T_ast <= 100)
                    muG=40.785*Fc*((MW*T)^(1/2))/(Omega_v*Vc^(2/3)); % Micropoise
                    value=muG*10^(-7); % Pa.s
                    % disp('### AVISO: Método de estimación en curso ###')
                else
                    disp('### AVISO: Temperatura fuera de rango ###')
                end
            end

        elseif strcmpi(propiedad,'Hf')
                value=-277.0*1000/(MW/1000);

            % Calor específico a presión constante, 'cp' [J/kg.K]
            % Perrys... no está en NASA
        elseif strcmpi(propiedad,'cp')
            C1=0.7454E+05; C2=2.5907E+05; C3=1.6073E+03; C4=1.732E+05; C5=712.4;
            if (200 <= T)&&(T <= 1500)
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
            else % Compruebo extrapolacion razonable...
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
                if muestra_aviso==1
                    disp('### [Butanol] AVISO: T fuera de rango en cálculo Cp_v. Se extrapolará ###')
                end  
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Conductividad vapor fuel (W/mK): Kee-CHEMKIN
            value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
              
        end
    end
    
    
    
    % #########################################################################
    % ______Propiedades ACETONA: Perry's Chem. Eng. Handbook, 8th Edition._____
    % #########################################################################
elseif strcmpi(comb,'acetona')   
    
    w= 0.307; % Factor acéntrico
    dipole= 2.9;  % Debye, [D]
    Vc= 209;    % [cm3/mol]
    
    % Propiedades ACETONA(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        % Presión de vapor del líquido, 'pvap' [Pa]
        if strcmpi(propiedad,'pvap')
            C1=69.006; C2=-5599.6; C3=-7.0985; C4=6.2237E-06; C5=2;
            if (178.45<=T)&&(T<=508.2)
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
                % Poco probable que T<178 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
        elseif strcmpi(propiedad,'Lv')
            C1=4.215E+07; C2=0.3397; C3=0; C4=0; C5=0;
            if (178.45 <= T)&&(T <= 508.2)
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
            else
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
                % Poco probable que T<178 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Densidad del líqido, 'rho' [kg/m3]
        elseif strcmpi(propiedad,'rho')
            C1=1.2332; C2=0.25886; C3=508.2; C4=0.2913;
            if (178.45 <= T)&&(T <= 508.2)
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
            else
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
                % Poco probable que T<178 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Viscosidad, 'mu' [Pa.s]
        elseif strcmpi(propiedad,'mu')
            C1=-14.918; C2=1023.4; C3=0.5961; C4=0; C5=0;
            if (190 <= T)&&(T <= 329.7)
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
                % Poco probable que T<190 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=135600; C2=-177; C3=0.2837; C4=0.000689; C5=0;
            if (178.45 <= T)&&(T <= 329.7)
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
            else
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
                % Poco probable que T<178 K...
                if muestra_aviso==1
                    disp('### AVISO: Se extrapolará ###')
                end
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Perrys Handbook. Relacion lineal con T, por lo que no parece
            % peligroso extrapolar...
            C1=0.2878; C2=-0.000427; C3=0; C4=0; C5=0;
            value=C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4;

        % Entalpía de formación, 'Hf' [J/Kg]
        % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
                value=-217.15*1000/(MW/1000);
        end
        
        
        % Propiedades ACETONA(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        if strcmpi(propiedad,'mu')
            C1=3.1005E-08; C2=0.9762; C3=23.139; C4=0;
            if (178.45 <= T)&&(T <= 1000)
                value=C1*(T^C2)/(1+C3/T+C4/(T^2));
            elseif T>1000 % Method of Chung et al., Poling et al., The Properties of Gases
                % and Liquids, pp. 9.7
                reduced_dipole=131.3*dipole/((Vc*Tc)^(1/2));
                Fc=1-0.2756*w+0.059035*reduced_dipole^4;
                T_ast=1.2593*Tr;
                A=1.16145; B=0.14874; C=0.52487; D=0.77320; E=2.16178; F=2.43787;
                Omega_v=(A*T_ast^(-B))+C*exp(-D*T_ast)+E*exp(-F*T_ast);
                if (0.3 <= T_ast)&&(T_ast <= 100)
                    muG=40.785*Fc*((MW*T)^(1/2))/(Omega_v*Vc^(2/3)); % Micropoise
                    value=muG*10^(-7); % Pa.s
                    % disp('### AVISO: Método de estimación en curso ###')
                else
                    disp('### AVISO: Temperatura fuera de rango ###')
                end
            end

        elseif strcmpi(propiedad,'Hf')
            value=-217.15*1000/(MW/1000);

            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=0.5704E+05; C2=1.632E+05; C3=1.6073E+03; C4=0.968E+05; C5=731.5;
            if (200 <= T)&&(T <= 1500)
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
            else % Compruebo extrapolacion razonable...
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
                if muestra_aviso==1
                    disp('### [Acetona] AVISO: T fuera de rango en cálculo Cp_v. Se extrapolará ###')
                end  
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Conductividad vapor fuel (W/mK): Kee-CHEMKIN
            value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
        end
    end
    
    
    
    % #########################################################################
    % _____Propiedades NAFTALENO: Perry's Chem. Eng. Handbook, 8th Edition.____
    % #########################################################################
elseif strcmpi(comb,'naftaleno')  
    
    w= 0.304; % Factor acéntrico
    dipole= 0;  % Debye, [D]
    Vc= 407;    % [cm3/mol]
    
    % Propiedades NAFTALENO(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        % Presión de vapor del líquido, 'pvap' [Pa]
        if strcmpi(propiedad,'pvap')
            C1=62.964; C2=-8137.5; C3=-5.6317; C4=2.2675E-18; C5=6;
            if (353.43<=T)&&(T<=748.4)
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else % Compruebo que la extrapolación es razonable...
                value = exp(C1+C2/T+C3*log(T)+C4*T^C5);
                if muestra_aviso==1
                    disp('### [Naftaleno] AVISO: T fuera de rango en cálculo Pvap (< 353 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
        elseif strcmpi(propiedad,'Lv')
            C1=7.0911E+07; C2=0.46468; C3=0; C4=0; C5=0;
            if (353.43 <= T)&&(T <= 748.4)
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
            else % Compruebo que la extrapolación es razonable...
                value = (C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2+C5*Tr^3))/MW;
                if muestra_aviso==1
                    disp('### [Naftaleno] AVISO: T fuera de rango en cálculo Lv (< 353 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Densidad del líqido, 'rho' [kg/m3]
        elseif strcmpi(propiedad,'rho')
            C1=0.6348; C2=0.25838; C3=748.4; C4=0.27727;
            if (353.43 <= T)&&(T <= 748.4)
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
            else % Compruebo que la extrapolación es razonable...
                value=(C1/C2^(1+(1-T/C3)^C4))*MW;
                if muestra_aviso==1
                    disp('### [Naftaleno] AVISO: T fuera de rango en cálculo rho_l (< 353 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Viscosidad, 'mu' [Pa.s]
        elseif strcmpi(propiedad,'mu')
            C1=-19.308; C2=1822.5; C3=1.218; C4=0; C5=0;
            if (353.43 <= T)&&(T <= 633.15)
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
            else % Compruebo que la extrapolación es razonable...
                value=exp(C1+C2/T+C3*log(T)+C4*T^C5);
                if muestra_aviso==1
                    disp('### [Naftaleno] AVISO: T fuera de rango en cálculo mu_l (< 353 K = T_fus). Se extrapolará ###')
                end
                
            end
            
            % Calor específico liq., 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=29800; C2=527.5; C3=0; C4=0; C5=0;
            if (353.43 <= T)&&(T <= 491.14)
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
            else  % Compruebo que la extrapolación es razonable...
                value=(C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4)/MW;
                if muestra_aviso==1
                    disp('### [Naftaleno] AVISO: T fuera de rango en cálculo C_l (< 353 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Perrys Handbook. Relacion lineal con T, por lo que no parece
            % peligroso extrapolar...
            C1=0.17096; C2=-0.00010059; C3=0; C4=0; C5=0;
            value=C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4;

        % Entalpía de formación, 'Hf' [J/Kg]
        % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
                value=150.58*1000/(MW/1000);
        end
        
        % Propiedades NAFTALENO(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        if strcmpi(propiedad,'mu')
            C1=6.4318E-07; C2=0.5389; C3=400.16; C4=0;
            if (353.43 <= T)&&(T <= 1000)
                value=C1*(T^C2)/(1+C3/T+C4/(T^2));
            elseif T>1000 % Method of Chung et al., Poling et al., The Properties of Gases
                % and Liquids, pp. 9.7
                reduced_dipole=131.3*dipole/((Vc*Tc)^(1/2));
                Fc=1-0.2756*w+0.059035*reduced_dipole^4;
                T_ast=1.2593*Tr;
                A=1.16145; B=0.14874; C=0.52487; D=0.77320; E=2.16178; F=2.43787;
                Omega_v=(A*T_ast^(-B))+C*exp(-D*T_ast)+E*exp(-F*T_ast);
                if (0.3 <= T_ast)&&(T_ast <= 100)
                    muG=40.785*Fc*((MW*T)^(1/2))/(Omega_v*Vc^(2/3)); % Micropoise
                    value=muG*10^(-7); % Pa.s
                    % disp('### AVISO: Método de estimación en curso ###')
                else
                    disp('### AVISO: Temperatura fuera de rango ###')
                end
            end

        elseif strcmpi(propiedad,'Hf')
                value=150.58*1000/(MW/1000);
                
            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            C1=0.6805E+05; C2=3.5494E+05; C3=1.4262E+03; C4=2.5984E+05; C5=650.1;
            if (200 <= T)&&(T <= 1500)
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
            else % Compruebo extrapolacion razonable...
                value=(C1 + C2*((C3/T)/sinh(C3/T))^2 + C4*((C5/T)/cosh(C5/T))^2)/MW;
                if muestra_aviso==1
                    disp('### [Naftaleno] AVISO: T fuera de rango en cálculo Cp_v. Se extrapolará ###')
                end  
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Conductividad vapor fuel (W/mK): Kee-CHEMKIN
            value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
        end
    end
    
    
% ############################################################################################
% Propiedades GLICERINA: HANDBOOK OF PHYSICAL PROPERTIES FOR HYDROCARBONS AND CHEMICALS (Yaws)
% ############################################################################################
elseif strcmpi(comb,'glicerina')
    
    w= 1.5937;       % Factor acéntrico (calculado) w=-log10[(Pvap@T=0.7*Tc)/Pc]-1
    dipole= 2.7;       % Debye, [D]
    Vc= 131.0435;    % [cm3/mol]

    
    % Propiedades GLICERINA(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        % Presión de vapor del líquido, 'pvap' [Pa]
        if strcmpi(propiedad,'pvap')
            if (290<=T)&&(T<=723)
                value = 133.322365*10^(-62.7929 -3.6585E+03/T + 3.4249E+01*log10(T) - 5.1940E-02*T + 2.2830E-05*T^2);
            else % Compruebo que la extrapolación es razonable...
                value = 133.322365*10^(-62.7929 -3.6585E+03/T + 3.4249E+01*log10(T) - 5.1940E-02*T + 2.2830E-05*T^2);
                if muestra_aviso==1
                    disp('### [Glicerina] AVISO: T fuera de rango en cálculo Pvap (< 290 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
        elseif strcmpi(propiedad,'Lv')
            if (290 <= T)&&(T <= 723)
                value = ((104.153*(1-T/723.00)^0.301)/MW)*1E+06;
            else % Compruebo que la extrapolación es razonable...
                value = ((104.153*(1-T/723.00)^0.301)/MW)*1E+06;
                if muestra_aviso==1
                    disp('### [Glicerina] AVISO: T fuera de rango en cálculo Lv (< 290 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Densidad del líqido, 'rho' [kg/m3]
        elseif strcmpi(propiedad,'rho')
            if (290 <= T)&&(T <= 723)
                value=1000*0.34908*0.24902^(-(1-T/723)^0.15410);
            else % Compruebo que la extrapolación es razonable...
                value=1000*0.34908*0.24902^(-(1-T/723)^0.15410);
                if muestra_aviso==1
                    disp('### [Glicerina] AVISO: T fuera de rango en cálculo rho_l (< 290 K = T_fus). Se extrapolará ###')
                end
            end

            % Viscosidad, 'mu' [Pa.s]
        elseif strcmpi(propiedad,'mu')
            A=-18.2152; B=4.2305E+03; C=2.8705E-02; D=-1.8648E-05;
            if (290 <= T)&&(T <= 723)
                value=1E-07*(10^(A+B/T+C*T+D*T^2));
            else % Compruebo que la extrapolación es razonable...
                value=1E-07*(10^(A+B/T+C*T+D*T^2));
                if muestra_aviso==1
                    disp('### [Glicerina] AVISO: T fuera de rango en cálculo mu_l (< 290 K = T_fus). Se extrapolará ###')
                end
            end
            
            % Calor específico liq., 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            if (290 <= T)&&(T <= 651)
                value=(132.145 + 8.6007E-01*T -1.9745E-03*T^2 +1.8068E-06*T^3)*(1000/MW);
            else % Compruebo que la extrapolación es razonable...
                value=(132.145 + 8.6007E-01*T -1.9745E-03*T^2 +1.8068E-06*T^3)*(1000/MW);
                if muestra_aviso==1
                    disp('### [Glicerina] AVISO: T fuera de rango en cálculo C_l (< 290 K = T_fus). Se extrapolará ###')
                end
            end

            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Relacion no-lineal (como en Perrys) pero casi... parece razonable extrapolar
            % cerca de los limites de Tmin y Tmax...
            value=10^(-0.3550 - 0.2097*(1-T/723)^(2/7)); % 293-550 K

        % Entalpía de formación, 'Hf' [J/Kg]
        % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
                value=-577.9*1000/(MW/1000);
        end
   
        % Propiedades GLICERINA(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        if strcmpi(propiedad,'mu')
            if (563 <= T)&&(T <= 993)
                value=(10^-7)*(-23.119 + 2.8879E+01*T - 3.4277E-05*T^2);
            else % Method of Chung et al., Poling et al., The Properties of Gases
                          % and Liquids, pp. 9.7
                reduced_dipole=131.3*dipole/((Vc*Tc)^(1/2));
                kappa=0.0682+4.704*3/MW;
                Fc=1-0.2756*w+0.059035*reduced_dipole^4+kappa;
                T_ast=1.2593*Tr;
                A=1.16145; B=0.14874; C=0.52487; D=0.77320; E=2.16178; F=2.43787;
                Omega_v=(A*T_ast^(-B))+C*exp(-D*T_ast)+E*exp(-F*T_ast);
                if (0.3 <= T_ast)&&(T_ast <= 100)
                    value=(40.785*Fc*((MW*T)^(1/2))/(Omega_v*Vc^(2/3)))*10^(-7); % Pa.s; 
                    % disp('### AVISO: Método de estimación en curso ###')
                else
                    disp('### AVISO: Temperatura fuera de rango ###')
                end
            end

        elseif strcmpi(propiedad,'Hf')
                value=-577.9*1000/(MW/1000);

            % Calor específico a presión constante, 'cp' [J/kg.K]
        elseif strcmpi(propiedad,'cp')
            if (298 <= T)&&(T <= 1200)
                value=(9.656 + 4.2826E-01*T - 2.6797E-04*T^2 + 3.1794E-08*T^3 + 2.7745E-11*T^4)*(1000/MW);
            else
                T_aux=1200;
                value=(9.656 + 4.2826E-01*T_aux - 2.6797E-04*T_aux^2 + 3.1794E-08*T_aux^3 + 2.7745E-11*T_aux^4)*(1000/MW);

                % OJO! Extrapolación NO parece razonable, Cp  se dispara a
                % altas temperaturas (>1200 K). Si tenemos que simular gicerina bajo estas 
                % condiciones --> Buscar nueva correlacion para Cp! (no aparece en Perrys, NIST...)
                % disp('### [Glicerina] AVISO: T fuera de rango en cálculo Cp_v (> 1200 K). NO se recomienda extrapolar por la función empleada ###')
                % disp('El programa se detendrá. Buscar nueva correlación para Cp_v de la glicerina ###')
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
             % Conductividad vapor fuel (W/mK): Kee-CHEMKIN
            value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
        end
    end
    
    
    % #########################################################################
    % ______________________Propiedades 1-METILNAFTALENO_______________________
    % #########################################################################
elseif strcmpi(comb,'1-metilnaftaleno')   
%     El caso del 1-MNP es especial, no se puede tratar como el
%   resto de combustibles al no disponerse de expresiones tan directas para
%   el cálculo de sus propiedades.
    
        % Propiedades no definidas en objeto gota que son necesarias para
    % cálculo de prop. del 1MNP:
    w= 0.348;     % Factor acéntrico, [adimensional]
    dipole= 0.5;  % Debye, [D]
    Vref= 139.37; % [cm3/mol]
    Tref= 293.15;
    Vc= 462;      % [cm3/mol]
    
    
    
    % Propiedades 1-METILNAFTALENO(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        % Presión de vapor del líquido, 'pvap' [Pa]
        % Poling et al., The Properties of Gases and Liquids (Ec. Antoine,
        % pp. 7.4; Parámetros ec. Antoine, pp. A.58)
        if strcmpi(propiedad,'pvap')
            A=4.16082; B=1826.948; C=195.002;
            if (389.93<=T)&&(T<=517.8) 
                value = (10^5)*(10^(A-B/(T+C-273.15)));
            else % Compruebo que la extrapolación es razonable...
                value = (10^5)*(10^(A-B/(T+C-273.15)));
                if muestra_aviso==1
                    disp('### [1-MNP] AVISO: T fuera de rango en cálculo Pvap (< 390 K). Se extrapolará ###')
                end
            end
            
            % Calor latente de vaporización del líquido, 'Lv' [J/kg]
            % Pitzer et al. corresponding-states correlation, Poling et al.,
            % The Properties of Gases and Liquids, pp. 7.18
        elseif strcmpi(propiedad,'Lv')
            tau=1-Tr;
            if ((0.6*Tc) <= T)&&(T <= (1*Tc))
                value = R*Tc*(7.08*tau^0.354 + 10.95*w*tau^0.456)*(1000/MW);
            else % Compruebo que la extrapolación es razonable...
                value = R*Tc*(7.08*tau^0.354 + 10.95*w*tau^0.456)*(1000/MW);
                if muestra_aviso==1
                    disp('### [1-MNP] AVISO: T fuera de rango en cálculo Lv (< 0.6*Tc=463 K). Se extrapolará ###')
                end
            end
            
            % Densidad del líqido, 'rho' [kg/m3]
            % Rackett, Yamada and Gunn eq., Poling et al., The Properties of
            % Gases and Liquids, pp. 4.35
        elseif strcmpi(propiedad,'rho')
            % No se indican límites de aplicabilidad de esta expresion,
            % pero parece bastante lineal en toda la zona de interés (250 K - Tb)
            phi=(1-T/Tc)^(2/7)-(1-Tref/Tc)^(2/7);
            V=Vref*(0.29056-0.08775*w)^phi;
            value=(MW/V)*1000;
            
            % Viscosidad, 'mu' [Pa.s]
            % Chemical Engineering Design, "Appendix C"
        elseif strcmpi(propiedad,'mu')
            VISA=862.89; VISB=361.76;
            value=(10^(VISA*(1/T-1/VISB)))/1000;
            
            
            % Calor específico a presión constante, 'cp' [J/kg.K]
            % Corresponding-states method, Bondi, refitted, Poling et al, The
            % Properties of Gases and Liquids, pp. 6.22
        elseif strcmpi(propiedad,'cp')
            if T<=(0.99*Tc)
                % Calculamos en primer lugar el calor específico a presión
                % constante del vapor (lo necesitamos en J/mol.K)
                a0=-5.637; a1=98.625E-03; a2=-4.956E-05; a3=-1.033E-08; a4=1.281E-11;
                if (200 <= T)&&(T <= 1500)
                    CpG=R*(a0+a1*T+a2*T^2+a3*T^3+a4*T^4);
                else
                    disp('### AVISO: Temperatura fuera de rango ###')
                end
                value=(CpG+R*(1.586+0.49/(1-Tr)+w*(4.2775+(6.3/Tr)*(1-Tr)^(1/3)+0.4355/(1-Tr))))*(1000/MW);
            else
                % A efectos prácticos cubre todo el rango alcanzable en el liquido...
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
            % Latini et al. method, Poling et al., The Properties of Gases and
            % Liquids, pp. 10.44
        elseif strcmpi(propiedad,'lambda')
            % Válido para T<=Tb (se da por hecho de que es así...). Además,
            % se comprueba que la dependencia es bastante lineal
            
            % Para aromatics:
            A_ast=0.0346; alpha=1.2; betta=1; gamma=0.167;
            % Por lo tanto:
            A=A_ast*(Tb^alpha)/((MW^betta)*(Tc^gamma));
            value=(A*(1-Tr)^0.38)/(Tr^(1/6));
        
        % Entalpía de formación, 'Hf' [J/Kg]
        % Debe ser el mismo que el 'vapor' por cómo se calcular h_i
        elseif strcmpi(propiedad,'Hf')
                value=116.9*1000/(MW/1000);
        end
        
        % Propiedades 1-METILNAFTALENO(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        
        % Viscosidad, 'mu' [Pa.s]
        % Method of Chung et al., Poling et al., The Properties of Gases
        % and Liquids, pp. 9.7
        if strcmpi(propiedad,'mu')
            reduced_dipole=131.3*dipole/((Vc*Tc)^(1/2));
            Fc=1-0.2756*w+0.059035*reduced_dipole^4;
            T_ast=1.2593*Tr;
            A=1.16145; B=0.14874; C=0.52487; D=0.77320; E=2.16178; F=2.43787;
            Omega_v=(A*T_ast^(-B))+C*exp(-D*T_ast)+E*exp(-F*T_ast);
            if (0.3 <= T_ast)&&(T_ast <= 100)
                muG=40.785*Fc*((MW*T)^(1/2))/(Omega_v*Vc^(2/3)); % Micropoise
                value=muG*10^(-7); % Pa.s
            else
                disp('### AVISO: Temperatura fuera de rango ###')
            end
            
        elseif strcmpi(propiedad,'Hf')
                value=116.9*1000/(MW/1000);

            % Calor específico a presión constante, 'cp' [J/kg.K]
            % Poling et al., The Properties of Gases and Liquids, pp. A.44
        elseif strcmpi(propiedad,'cp')
            a0=-5.637; a1=98.625E-03; a2=-4.956E-05; a3=-1.033E-08; a4=1.281E-11;
            if (273 <= T)&&(T <= 1200) 
                value=R*(a0+a1*T+a2*T^2+a3*T^3+a4*T^4)*(1000/MW);
            else  % OJO! Extrapolación NO parece razonable, Cp  se dispara a
                % altas temperaturas (>1200 K). Si tenemos que simular gicerina bajo estas 
                % condiciones --> Buscar nueva correlacion para Cp! (no aparece en Perrys, NIST...)
                disp('### [1-MNP] AVISO: T fuera de rango en cálculo Cp_v (> 1200 K). NO se recomienda extrapolar por la función empleada ###')
                disp('El programa se detendrá. Buscar nueva correlación para Cp_v del 1-MNP ###')
            end
            
            % Conductividad térmica, 'lambda' [W/m.K]
        elseif strcmpi(propiedad,'lambda')
            % Conductividad vapor fuel (W/mK): Kee-CHEMKIN
            value= calcula_lambda_vapor_Kee(comb, obj, T); %(* W/m/K *)
            
        end
    end
end



end

