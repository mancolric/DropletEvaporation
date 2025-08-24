function lambda = calcula_lambda_vapor_Kee(compuesto, objeto, T)
% Función que devuelve lambda vapor fuel/gas con método Kee (teria cinet. gases)
% Con los parámetros entrantes, puede calcular la conductividad tanto de:
% - fuel monocomponente: 'compuesto' es heptano/octano/etanol, etc.
% - gas monocomponente: 'compuesto' es N2/CO2, etc.
% Importante hacer esta distinción porque el método de cálculo me requiere
% calcular el calor específico a volumen cte (Cv) del compuesto en cuestión

R=8.314; %(J/molK)
% Listado de combustibles y gases que tenemos, de forma que podamos saber
% qué estamos calculando:
combustibles = {'heptano' 'hexano' 'octano' 'dodecano' 'tetradecano' 'hexadecano' 'eicosano' 'metanol' 'etanol' 'butanol' 'acetona' 'naftaleno' '1-metilnaftaleno' 'glicerina'};
gases= {'N2' 'H2O' 'CO2' 'O2'};

ind_f=find(strcmpi(combustibles,compuesto), 1);
ind_g=find(strcmpi(gases,compuesto), 1);

if ~isempty(ind_f) % Estoy calculando lambda de un fuel
    
    % OJO! Puede que tenga varios combustibles... debo sacar el indice del
    % fuel actual:
    for iD_combustible=1:size(objeto.tipo_combustible,2)
        if strcmpi(objeto.tipo_combustible{iD_combustible}, compuesto)
            indice_comb=iD_combustible;
            break
        end
    end
    
    
    % El "objeto" entrante es el objeto gota
    mw=objeto.mw(indice_comb)*1000; %(kg/kmol)
    sigma=objeto.sigma(indice_comb); %(m)
    e_k=objeto.e_k(indice_comb); %(K)
    Zrot_298=objeto.Zrot_298(indice_comb); %(-)
    molecula_lineal=objeto.molecula_lineal(indice_comb);
    
    % Puedo calcular el Cp del vapor de combustible a la temperatura T con:
    Cp=devuelve_propiedad_fuel(objeto, 'cp', compuesto, 'vapor', T); %J/kg/k

    
elseif ~isempty(ind_g) % Estoy calculando lambda de un gas
    % El "objeto" entrante no sirve de nada (viene vacío). Obtengo las
    % propiedades cinéticas del gas de la funcion
    % propiedades_Tcinetica_gasMonocomp:
    propiedades_Tcinetica=propiedades_Tcinetica_gasMonocomp(compuesto);
    mw=propiedades_Tcinetica(1)*1000; %(kg/kmol)
    sigma=propiedades_Tcinetica(2); %(m)
    e_k=propiedades_Tcinetica(3); %(K)
    Zrot_298=propiedades_Tcinetica(4); %(-)
    molecula_lineal=propiedades_Tcinetica(5); %(-)
    
    % Puedo calcular el Cp del vapor de gas monocomp. a la temperatura T con:
    Cp=calcula_Cp_gas_monocomponente(compuesto, T); %J/kg/k 
    
else % No sé qué estoy calculando. Avisar error
    disp('ERROR: NO SE ENCUENTRA EL COMPUESTO (cálculo conductividad vapor por teoría cinética)')
end

% Primero de todo calculo el Cv a partir de Cp:
Cp_mol = Cp*mw/1000; %J/molK
% Y el Cv sería (gas ideal: Cp=Cv+R):
Cv = (Cp_mol - R); %(J/molK)


% Correlacion de Kee (2):
k_Boltz=1.38065e-23;
num_Avog=6.022140857e23;
m_reduc=mw/1000/num_Avog; %(kg)
T_est=T/e_k;
a1=1.0548; a2=0.15504; a3=0.55909; a4=2.1705;
omega_1 = a1*T_est^(-a2) +(T_est+a3)^(-a4);
b1=1.0413; b2=0.11930; b3=0.43628; b4=1.6041;
omega_2 = b1*T_est^(-b2) +(T_est+b3)^(-b4);

mu_k=5/16*sqrt(pi*m_reduc*k_Boltz*T)/(pi*sigma^2*omega_2);

D_kk=3/8*sqrt(pi*k_Boltz^3*T^3/m_reduc)/(101325*pi*sigma^2*omega_1);

if molecula_lineal==1
    Cv_trans=3/2*R;
    Cv_rot=R;
    Cv_vib=Cv-5/2*R;
else
    Cv_trans=3/2*R;
    Cv_rot=3/2*R;
    Cv_vib=Cv-3*R;
end

F_1=1 + 0.5*pi^1.5*(e_k/T)^0.5 +(0.25*pi^2+2)*(e_k/T) + pi^1.5*(e_k/T)^1.5;
F_2=1 + 0.5*pi^1.5*(e_k/298)^0.5 +(0.25*pi^2+2)*(e_k/298) + pi^1.5*(e_k/298)^1.5;
Z_rotac = Zrot_298*F_2/F_1;

rho=101325*(mw/1000)/(R*T);

A=5/2*rho*D_kk/mu_k;
B=Z_rotac + 2/pi*(5/3*Cv_rot/R + rho*D_kk/mu_k);

f_trans=5/2*(1-(2*Cv_rot*A)/(pi*Cv_trans*B));
f_rot=rho*D_kk/mu_k*(1+2*A/(pi*B));
f_vib=rho*D_kk/mu_k;

lambda=mu_k/(mw/1000)*(f_trans*Cv_trans + f_rot*Cv_rot + f_vib*Cv_vib);

if (strcmpi('H2O',compuesto))
  lambda=(2.0103 - 7.9139*T/1000 + 35.922*(T/1000)^2 - 41.39*(T/1000)^3 +  35.993*(T/1000)^4 - 18.974*(T/1000)^5 + 4.1531*(T/1000)^6)*1/10^2; %Perrys (* W/m/K *)
    
end

return

