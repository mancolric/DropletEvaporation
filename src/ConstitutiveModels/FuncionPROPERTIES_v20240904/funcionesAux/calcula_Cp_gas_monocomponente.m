function Cp_g = calcula_Cp_gas_monocomponente(tipo_gas, T)
%Devuelve el Cp del gas especificado a la temperatura T
R=8.314;

if strcmpi(tipo_gas,'N2')
    mw_g = 28; %(* kg/kmol *)
    % NASA polynomial fit (J/kgK)
    if (200 <= T)&&(T <= 1000) %Primer intervalo
        a1=3.53100528E+0; a2=-1.23660987E-4; a3=-5.02999437E-07; a4=2.43530612E-09; a5=-1.40881235E-12;
        Cp_g=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/mw_g;
    elseif (T > 1000) %Segundo intervalo
        a1=2.95257626E+00; a2=1.39690057E-3; a3=-4.92631691E-07; a4=7.86010367E-11; a5=-4.60755321e-15;
        Cp_g=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/mw_g;
    end
    
elseif strcmpi(tipo_gas,'O2')
    mw_g = 32; %(* kg/kmol *)
    % NASA polynomial fit (J/kgK)
    if (200 <= T)&&(T <= 1000) %Primer intervalo
        a1=3.78245636e0; a2=-2.99673415e-3; a3=9.84730200e-6; a4=-9.68129508e-9; a5=3.24372836e-12;
        Cp_g=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/mw_g;
    elseif (T > 1000) %Segundo intervalo
        a1=3.66096083e0; a2=6.56365523e-4; a3=-1.41149485e-7; a4=2.05797658e-11; a5=-1.29913248e-15;
        Cp_g=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/mw_g;
    end
    
elseif strcmpi(tipo_gas,'CO2')
    mw_g = 44;  %(* kg/kmol *)
    if (200 <= T)&&(T <= 1000) %Primer intervalo
        a1=2.35677352e0; a2=8.98459677e-3; a3=-7.12356269e-6; a4=2.45919022e-9; a5=-1.43699548e-13;
        Cp_g=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/mw_g;
    elseif (T > 1000) %Segundo intervalo
        a1=4.63659493e0; a2=2.74131991e-3; a3=-9.95828531e-7; a4=1.60373011e-10; a5=-9.16103468e-15;
        Cp_g=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/mw_g;
    end
    
elseif strcmpi(tipo_gas,'H2O')
    mw_g = 18;  %(* kg/kmol *)
    if (200 <= T)&&(T <= 1000) %Primer intervalo
        a1=4.19864056e0; a2=-2.03643410e-3; a3=6.52040211e-6; a4=-5.48797062e-9; a5=1.77197817e-12;
        Cp_g=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/mw_g;
    elseif (T > 1000) %Segundo intervalo
        a1=2.67703787e0; a2=2.97318329e-3; a3=-7.73769690e-7; a4=9.44336689e-11; a5=-4.26900959e-15;
        Cp_g=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/mw_g;
    end


elseif strcmpi(tipo_gas,'CO')
    % http://combustion.berkeley.edu/gri-mech/version30/files30/thermo30.dat
    % (NASA)
    mw_g = 28;  %(* kg/kmol *)
    if (200 <= T)&&(T <= 1000) %Primer intervalo
        a1= 3.57953347E+00; a2= -6.10353680E-04; a3= 1.01681433E-06; a4= 9.07005884E-10; a5= -9.04424499E-13;
        Cp_g=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/mw_g;
    elseif (T > 1000) %Segundo intervalo
        a1= 2.71518561E+00; a2= 2.06252743E-03; a3= -9.98825771E-07; a4= 2.30053008E-10; a5= -2.03647716E-14;
        Cp_g=R*1000*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4)/mw_g;
    end
end

end

