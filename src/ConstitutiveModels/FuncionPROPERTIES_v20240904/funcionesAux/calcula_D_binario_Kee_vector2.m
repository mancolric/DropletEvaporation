function D_ab = calcula_D_binario_Kee_vector2(Compuesto1, Compuesto2, T_eval, epsilon_f_eval)

if isobject(Compuesto1)
    N_puntos    = size(T_eval,2);
    D_ab        = zeros(1, N_puntos);

    prop_gas    = propiedades_Tcinetica_gasMonocomp(Compuesto2); %[mw, sigma, e_k, Zrot_298, molecula_lineal]
    mw_g        = prop_gas(1)*1000; %(kg/kmol)
    sigma_g     = prop_gas(2); %(m)
    ek_g        = prop_gas(3); %(K)
    
    mw_f        = Compuesto1.mw*1000*epsilon_f_eval;
    sigma_f     = Compuesto1.sigma*epsilon_f_eval;
    ek_f        = Compuesto1.e_k*epsilon_f_eval;

    k_Boltz     = 1.38065e-23;
    num_Avog    = 6.022140857e23;
    M_ab        = mw_f.*mw_g./(mw_f+mw_g);
    m_reduc_f   = mw_f./1000/num_Avog; %(kg)
    m_reduc_g   = mw_g./1000/num_Avog; %(kg)
    m_reduc_fg  = (m_reduc_g.*m_reduc_f)./(m_reduc_g+m_reduc_f);
    sigma_ab    = 0.5.*(sigma_g+sigma_f); %(m)
    ek_fg       = (ek_f.*ek_g).^0.5;
    T_est       = T_eval./ek_fg;
    a1=1.0548; a2=0.15504; a3=0.55909; a4=2.1705;
    omega_d     = a1.*T_est.^(-a2) +(T_est+a3).^(-a4);

    D_ab        = 3/16*((2.*pi.*k_Boltz.^3.*T_eval.^3./m_reduc_fg).^0.5)./(101325*pi.*sigma_ab.^2.*omega_d);

elseif ischar(Compuesto1)

    N_puntos    = size(T_eval,2);
    D_ab        = zeros(1, N_puntos);

    prop_gas2   = propiedades_Tcinetica_gasMonocomp(Compuesto2); %[mw, sigma, e_k, Zrot_298, molecula_lineal]
    mw_g        = prop_gas2(1)*1000; %(kg/kmol)
    sigma_g     = prop_gas2(2); %(m)
    ek_g        = prop_gas2(3); %(K)

    prop_gas1   = propiedades_Tcinetica_gasMonocomp(Compuesto1); %[mw, sigma, e_k, Zrot_298, molecula_lineal]
    mw_f        = prop_gas1(1)*1000; %(kg/kmol)
    sigma_f     = prop_gas1(2); %(m)
    ek_f        = prop_gas1(3); %(K)

    k_Boltz     = 1.38065e-23;
    num_Avog    = 6.022140857e23;
    M_ab        = mw_f.*mw_g./(mw_f+mw_g);
    m_reduc_f   = mw_f./1000/num_Avog; %(kg)
    m_reduc_g   = mw_g./1000/num_Avog; %(kg)
    m_reduc_fg  = (m_reduc_g.*m_reduc_f)./(m_reduc_g+m_reduc_f);
    sigma_ab    = 0.5.*(sigma_g+sigma_f); %(m)
    ek_fg       = (ek_f.*ek_g).^0.5;
    T_est       = T_eval./ek_fg;
    a1=1.0548; a2=0.15504; a3=0.55909; a4=2.1705;
    omega_d     = a1.*T_est.^(-a2) +(T_est+a3).^(-a4);

    D_ab        = 3/16*((2.*pi.*k_Boltz.^3.*T_eval.^3./m_reduc_fg).^0.5)./(101325*pi.*sigma_ab.^2.*omega_d);

else

    error('Compuesto 1 en calculo de Difusividad binaria no es ni inerte ni fuel o Fuel metido en Compuesto2')
end
return