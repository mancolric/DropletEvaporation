classdef clase_gota<handle
    
    properties
        
        % PROPIEDADES INSTRÍNSECAS DEL COMBUSTIBLE/S: SE DEFINEN EN EL BUILDER
        tipo_combustible % hexadecano, etanol, butanol, naftaleno...
        mw %Masa molecular (kg/mol)
        Tf %Temp fusión a Patm(K)
        Tb %Temp ebullicion a Patm (K)
        Tc %Temp critica (K)
        Pc %Presion critica (Pa)
        sigma % Lennard-Jones length para propiedades de trasnporte (m)
        e_k % %Lennard-Jones energy para propiedades de trasnporte (K)
        Zrot_298 %Z rotacional a 298K para prop. transporte (-). Para moléculas grandes (todos comb.) vale 1. Para moléculas pequeñas como N2 o H2O puede valer >1
        molecula_lineal % 1 si la molecula es lineal, 0 si no lo es (para cálculo lambda_v)
        
        
        % PROPIEDADES DEPENDIENTES DE T: EN actualiza_gota
        Pvap %Presión de vapor (Pa)
        rhol %Densidad liquido(kg/m3)
        Cpl %Calor especifico liquido (J/kgK)
        Cpv %Calor especifico vapor (J/kgK)
        Lv %Entalpía de vaporizacion (J/kg)
        mu_v %Viscosidad vapor (Pa·s)
        mu_l %Viscosidad líquido (Pa·s)
        lambda_v %Conductividad termica vapor (W/mK)
        lambda_l %Conductividad termica liquido (W/mK)
        Dif_l %Difusividad del líquido (m2/s)
        
        % PROPIEDADES DE EVALUACIÓN: Indica en qué condiciones se tiene
        % actualizado el objeto gota. TAMBIÉN EN actualiza_gota
        lugar_evaluacion % liquido, vapor...
        T_evaluacion %Temperatura (K)
        Y_evaluacion % Fracción másica de cada combustible (-)
    end
    
    methods
        
        %"Builder":
        function obj=clase_gota(combustible)
            % En el Builder no defino más que las propiedades que
            % van a mantenerse constantes durante toda la ejecución: MW, Tb,
            % Tc, etc. (constantes del combustible).
            
            for i=1:length(combustible)
                if strcmpi(combustible(i),'Octano') 
                    obj.tipo_combustible{i}='Octano';
                    obj.mw(i)=114.229/1000;
                    obj.Tf(i)=216.39;
                    obj.Tb(i)=398.82;
                    obj.Tc(i)=568.70;
                    obj.Pc(i)=2490000;
                    obj.sigma(i)=6.170*1e-10; %San Diego Mec
                    obj.e_k(i)=494.0; %San Diego Mec
                    obj.Zrot_298(i)=1; %San Diego Mec
                    obj.molecula_lineal(i)=1;
                    
                elseif strcmpi(combustible(i),'Heptano')
                    obj.tipo_combustible{i}='Heptano';
                    obj.mw(i)=100.2019/1000;
                    obj.Tf(i)=182.59;
                    obj.Tb(i)=371.57;
                    obj.Tc(i)=540.20;
                    obj.Pc(i)=2740000;
                    obj.sigma(i)=6.253*1e-10; %San Diego Mec
                    obj.e_k(i)=459.6; %San Diego Mec
                    obj.Zrot_298(i)=1; %San Diego Mec
                    obj.molecula_lineal(i)=1;
                    
                    elseif strcmpi(combustible(i),'Hexano')
                    obj.tipo_combustible{i}='Hexano';
                    obj.mw(i)=86.1754/1000;
                    obj.Tf(i)=178;
                    obj.Tb(i)=341.9;
                    obj.Tc(i)=507.6;
                    obj.Pc(i)=30.2e5;
                    obj.sigma(i)=5.946*1e-10; %San Diego Mec
                    obj.e_k(i)=427.4; %San Diego Mec
                    obj.Zrot_298(i)=1; %San Diego Mec
                    obj.molecula_lineal(i)=1;
                    
                elseif strcmpi(combustible(i),'Etanol')
                    obj.tipo_combustible{i}='Etanol';
                    obj.mw(i)=46.068/1000;
                    obj.Tf(i)=159.05;
                    obj.Tb(i)=351.80;
                    obj.Tc(i)=513.92;
                    obj.Pc(i)=6148000;
                    obj.sigma(i)=4.530*1e-10; %Poling (Appendix B) +artículo Sazhin + CHEMKIN
                    obj.e_k(i)=362.6; % %Poling (Appendix B) +artículo Sazhin + CHEMKIN
                    obj.Zrot_298(i)=1; %CHEMKIN
                    obj.molecula_lineal(i)=1;

                elseif strcmpi(combustible(i),'Metanol') 
                    obj.tipo_combustible{i}='Metanol';
                    obj.mw(i)=32.04/1000;
                    obj.Tf(i)=175.49;
                    obj.Tb(i)=337.69;
                    obj.Tc(i)=512.64;
                    obj.Pc(i)=8097000;
                    obj.sigma(i)=3.626*1e-10; %Poling (Appendix B) + San Diego Mech + CHEMKIN
                    obj.e_k(i)=481.8; % %Poling (Appendix B) + San Diego Mech + CHEMKIN
                    obj.Zrot_298(i)=1; % %CHEMKIN
                    obj.molecula_lineal(i)=0; %CHEMKIN
    
                    
                elseif strcmpi(combustible(i),'Butanol') 
                    obj.tipo_combustible{i}='Butanol';
                    obj.mw(i)=74.123/1000;
                    obj.Tf(i)=183.35;
                    obj.Tb(i)=390.6;
                    obj.Tc(i)=562.93;
                    obj.Pc(i)=45*1e5;
                    obj.sigma(i)=2.393*(562.93/(45*1e5/101325))^(1/3)*1e-10; % De thumb rule Kee (pag. 520 pdf)
                    obj.e_k(i)=1.18*390.6;  % De thumb rule Kee (pag. 520 pdf): e_k=1.18*Tb
                    obj.Zrot_298(i)=1;
                    obj.molecula_lineal(i)=1;

                
                elseif strcmpi(combustible(i),'Acetona')
                    obj.tipo_combustible{i}='Acetona';
                    obj.mw(i)=58.079/1000;
                    obj.Tf(i)=178.50;
                    obj.Tb(i)=329.65;
                    obj.Tc(i)=508.2;
                    obj.Pc(i)=4701000;
                    obj.sigma(i)=4.600*1e-10; % Poling (Appendix B)
                    obj.e_k(i)=560.2;  % Poling (Appendix B)
                    obj.Zrot_298(i)=1;
                    obj.molecula_lineal(i)=1;
                    
                elseif strcmpi(combustible(i),'Dodecano')
                    obj.tipo_combustible{i}='Dodecano';
                    obj.mw(i)=170.34/1000;
                    obj.Tf(i)=263.57;
                    obj.Tb(i)=489.5;
                    obj.Tc(i)=658.3;
                    obj.Pc(i)=1820000;
                    obj.sigma(i)=2.393*(658.3/(1820000/101325))^(1/3)*1e-10; % De thumb rule Kee (pag. 520 pdf)
                    obj.e_k(i)=1.18*489.5;  % De thumb rule Kee (pag. 520 pdf): e_k=1.18*Tb
                    obj.Zrot_298(i)=1; % Molécula más grande que heptano, que ya tiene Z=1..
                   	obj.molecula_lineal(i)=1;
                    
                    elseif strcmpi(combustible(i),'Tetradecano')
                    obj.tipo_combustible{i}='Tetradecano';
                    obj.mw(i)=198.388/1000;
                    obj.Tf(i)=278.7;
                    obj.Tb(i)=523;
                    obj.Tc(i)=693;
                    obj.Pc(i)=16e5;
                    obj.sigma(i)=2.393*(693/(16e5/101325))^(1/3)*1e-10; % De thumb rule Kee (pag. 520 pdf)
                    obj.e_k(i)=1.18*523;  % De thumb rule Kee (pag. 520 pdf): e_k=1.18*Tb
                    obj.Zrot_298(i)=1; % Molécula más grande que heptano, que ya tiene Z=1.

                elseif strcmpi(combustible(i),'Hexadecano') 
                    obj.tipo_combustible{i}='Hexadecano';
                    obj.mw(i)=226.448/1000;
                    obj.Tf(i)=291.32;
                    obj.Tb(i)=560;
                    obj.Tc(i)=717;
                    obj.Pc(i)=1420000;
                    obj.sigma(i)=2.393*(717/(1420000/101325))^(1/3)*1e-10; % De thumb rule Kee (pag. 520 pdf)
                    obj.e_k(i)=1.18*560;  % De thumb rule Kee (pag. 520 pdf)
                    obj.Zrot_298(i)=1; % Molécula más grande que heptano, que ya tiene Z=1..
                   	obj.molecula_lineal(i)=1;
                    
                elseif strcmpi(combustible(i),'Eicosano')
                    obj.tipo_combustible{i}='Eicosano';
                    obj.mw(i)=282.556/1000;
                    obj.Tf(i)=309.95;
                    obj.Tb(i)=617;
                    obj.Tc(i)=767;
                    obj.Pc(i)=1110000;
                    obj.sigma(i)=2.393*(767/(1110000/101325))^(1/3)*1e-10; % De thumb rule Kee (pag. 520 pdf)
                    obj.e_k(i)=1.18*617;  % De thumb rule Kee (pag. 520 pdf)
                    obj.Zrot_298(i)=1; % Molécula más grande que heptano, que ya tiene Z=1..
                   	obj.molecula_lineal(i)=1;
                    
                elseif strcmpi(combustible(i),'Naftaleno') 
                    obj.tipo_combustible{i}='Naftaleno';
                    obj.mw(i)=128.174/1000;
                    obj.Tf(i)=351.35;
                    obj.Tb(i)=491.1;
                    obj.Tc(i)=748.4;
                    obj.Pc(i)=4050000;
                    obj.sigma(i)=2.393*(748.4/(4050000/101325))^(1/3)*1e-10; % De thumb rule Kee (pag. 520 pdf)
                    obj.e_k(i)=1.18*491.1;  % De thumb rule Kee (pag. 520 pdf)
                    obj.Zrot_298(i)=1;
                   	obj.molecula_lineal(i)=1;
                    
                elseif strcmpi(combustible(i),'1-Metilnaftaleno') 
                    obj.tipo_combustible{i}='1-Metilnaftaleno';
                    obj.mw(i)=142.201/1000;
                    obj.Tf(i)=242.69;
                    obj.Tb(i)=517.8;
                    obj.Tc(i)=772;
                    obj.Pc(i)=36e5; %NIST
                    obj.sigma(i)=2.393*(772/(36e5/101325))^(1/3)*1e-10; % De thumb rule Kee (pag. 520 pdf)
                    obj.e_k(i)=1.18*517.8;  % De thumb rule Kee (pag. 520 pdf)
                    obj.Zrot_298(i)=1; 
                   	obj.molecula_lineal(i)=1;
                    
                    elseif strcmpi(combustible(i),'Glicerina')
                    obj.tipo_combustible{i}='Glicerina';
                    obj.mw(i)=92.0938/1000;
                    obj.Tf(i)=290;
                    obj.Tb(i)=560.3;
                    obj.Tc(i)=850;
                    obj.Pc(i)=75*1e5;
                    obj.sigma(i)=2.393*(850/(75*1e5/101325))^(1/3)*1e-10; % De thumb rule Kee (pag. 520 pdf)
                    obj.e_k(i)=1.18*560.3;  % De thumb rule Kee (pag. 520 pdf)
                    obj.Zrot_298(i)=1; 
                   	obj.molecula_lineal(i)=1;
                    
                else
                    % Combustible no reconocido:
                    fprintf('ERROR: NO SE RECONOCE LA SIGUIENTE SUSTANCIA: %s\n', combustible{i})
                    fprintf('Consultar lista de combustibles disponibles en función clase_gota\n')
                    fprintf('Se suspende la ejecución del programa \n');
                    stop %Para forzar a suspender ejecución del código
               
                end
            end
            
        end
        
                
        
        
        %##################################################################
        % ACTUALIZAR_GOTA:
        %Método para actualizar propiedades del comustible cuando cambia la temperatura o composición de la gota:
        function actualizar_gota(obj,vec_radio, nueva_temp,nueva_fracc,lugar_evaluacion)
        % OJO! Aquí es clave el lugar_evaluacion para saber qué está entrando y qué propiedades se van a actualizar:    

        % -'liquido_bulk' --> se actualizan rhol, Cpl, lambda_l y mu_l, todas
        % ellas calculadas a la T y composición integrales medias

        % -'liquido_superficie' --> se actualizan Lv y Pvap, calculadas a
        % la T y composición de la superficie (último valor de los arrays)
        
        % -'vapor' --> se actualizan Cpv, mu_v y lambda_v. Los valores de
        % nueva_T y nueva_Y son valores únicos y no perfiles radiales.
        
       if strcmpi(lugar_evaluacion, 'liquido_bulk')
            obj.lugar_evaluacion='liquido_bulk';
           
            % Calculo temperatura y composición integral promedio:
             T_eval=temperatura_integral_gota(obj, vec_radio,nueva_temp, nueva_fracc);
             y_eval=fraccion_integral_gota(obj, vec_radio,nueva_fracc, nueva_temp);
             
             % Actualizo tan solo las propiedades afectadas:
             obj.rhol=calcula_rhol_fuel(obj, y_eval,T_eval,lugar_evaluacion);
             obj.Cpl = calcula_Cp_fuel(obj,y_eval,T_eval,lugar_evaluacion); % calcula_Cp se encarga tanto de Cpl como de Cpv
             obj.lambda_l = calcula_lambda_fuel(obj,y_eval,T_eval,lugar_evaluacion); % calcula_lambda se encarga tanto de lambda_l como de lambda_v
             obj.mu_l=calcula_mu_fuel(obj,y_eval,T_eval,lugar_evaluacion); % calcula_mu se encarga tanto de mu_l como de mu_v
             obj.Dif_l=calcula_Dliq_fuel(obj, y_eval,T_eval,lugar_evaluacion);
             
             % Y elimino valores previos (de otra evaluación anterior) para
             % evitar confusiones:
             obj.Lv=[]; obj.Pvap=[]; obj.Cpv=[];obj.mu_v=[]; obj.lambda_v=[];
             

       elseif strcmpi(lugar_evaluacion, 'liquido_superficie')
            obj.lugar_evaluacion='liquido_superficie';

              %Obtengo temperatura y composición en la superficie:
            T_eval=nueva_temp(end); %Ts
            y_eval=nueva_fracc(:,end); %Ys  

            % Actualizo tan solo las propiedades afectadas:
            obj.Lv = calcula_Lv_fuel(obj, 0, T_eval, lugar_evaluacion); 
            obj.Pvap=calcula_Pvap_fuel(obj, 0, T_eval, lugar_evaluacion);
            
            % Y elimino valores previos (de otra evaluación anterior) para
             % evitar confusiones:
             obj.Cpv=[];obj.mu_v=[]; obj.lambda_v=[]; obj.rhol=[]; obj.Cpl=[]; obj.lambda_l=[]; obj.mu_l=[]; obj.Dif_l=[];
           
       elseif  strcmpi(lugar_evaluacion, 'vapor')
           obj.lugar_evaluacion='vapor';
           
           % En este caso no me llega perfil sino valores únicos de T y comp:
           T_eval=nueva_temp;
           y_eval=nueva_fracc;
           
           % Actualizo tan solo las propiedades afectadas:
           obj.Cpv =  calcula_Cp_fuel(obj,y_eval,T_eval,lugar_evaluacion); % calcula_Cp se encarga tanto de Cpl como de Cpv
           obj.mu_v = calcula_mu_fuel(obj,y_eval,T_eval,lugar_evaluacion); % calcula_mu se encarga tanto de mu_l como de mu_v
           obj.lambda_v = calcula_lambda_fuel(obj,y_eval,T_eval,lugar_evaluacion); % calcula_lambda se encarga tanto de lambda_l como de lambda_v
           
           % Y elimino valores previos (de otra evaluación anterior) para
             % evitar confusiones:
             obj.rhol=[]; obj.Cpl=[]; obj.lambda_l=[]; obj.mu_l=[]; obj.Lv=[]; obj.Pvap=[]; obj.Dif_l=[];
           
       else
           fprintf('ERROR EN LA EVALUACIÓN DE PROPIEDADES!\n')
           fprintf('Se suspende la ejecución del programa \n');
           stop %Para forzar a suspender ejecución del código
       end
                       
                  

       % Guardo (aunque es información redundante) la temperatura a
       % la que están evaluadas las propiedades actuales de la
       % gota. [Va bien para chequeos...]
       obj.T_evaluacion=T_eval;
       obj.Y_evaluacion=y_eval;
       
                
            
            
            
        end
        
    end
    
end