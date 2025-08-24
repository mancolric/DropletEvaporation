function D_ab = calcula_D_binario_Kee_vector(objeto, tipo_gas, T_eval, epsilon_f_eval)
% Función que devuelve D_ab de mezcla fuel-gas con metodo Kee-CHEMKIN
% Como D_ab es coeficiente binario, solo necesito saber:
% - Tipo de fuel (monoc)
% - Tipo de gas (monoc)
% - Temperatura

N_puntos=size(T_eval,2);
D_ab=zeros(1, N_puntos);

prop_gas= propiedades_Tcinetica_gasMonocomp(tipo_gas); %[mw, sigma, e_k, Zrot_298, molecula_lineal]
mw_g=prop_gas(1)*1000; %(kg/kmol)
sigma_g=prop_gas(2); %(m)
ek_g=prop_gas(3); %(K)
mw_f=objeto.mw*1000*epsilon_f_eval;
sigma_f=objeto.sigma*epsilon_f_eval;
ek_f=objeto.e_k*epsilon_f_eval;

k_Boltz=1.38065e-23;
num_Avog=6.022140857e23;
M_ab=mw_f.*mw_g./(mw_f+mw_g);
m_reduc_f=mw_f./1000/num_Avog; %(kg)
m_reduc_g=mw_g./1000/num_Avog; %(kg)
m_reduc_fg=(m_reduc_g.*m_reduc_f)./(m_reduc_g+m_reduc_f);
sigma_ab=0.5.*(sigma_g+sigma_f); %(m)
ek_fg=(ek_f.*ek_g).^0.5;
T_est=T_eval./ek_fg;
a1=1.0548; a2=0.15504; a3=0.55909; a4=2.1705;
omega_d = a1.*T_est.^(-a2) +(T_est+a3).^(-a4);

D_ab=3/16*((2.*pi.*k_Boltz.^3.*T_eval.^3./m_reduc_fg).^0.5)./(101325*pi.*sigma_ab.^2.*omega_d);

% for id_punto=1:N_puntos
% 1. Busco propiedades de teoría cinética del fuel.
% Si el fuel es multicomponente: ponderar por su fracción molar (Artículo
% Sazhin 2010, bicomponent...):
% mw_f=dot(epsilon_f_eval(:,id_punto), objeto.mw*1000); %(kg/kmol)
% sigma_f=dot(epsilon_f_eval(:,id_punto), objeto.sigma); %(m)
% ek_f=dot(epsilon_f_eval(:,id_punto), objeto.e_k); %(K)

% 2. Busco propiedades de teoría cinética de inertes:


% 3. Calculo coef. difusion binario entre ambas especies para el vector de temperaturas T_eval:   
% Correlacion de Kee:
% T=T_eval(id_punto);


 %(m2/s)

% end

return

