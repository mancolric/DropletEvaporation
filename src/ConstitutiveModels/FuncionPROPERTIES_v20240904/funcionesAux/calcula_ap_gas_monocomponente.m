function ap_g = calcula_ap_gas_monocomponente(tipo_gas, T)
%Devuelve el ap del gas especificado a la temperatura T

if strcmpi(tipo_gas,'N2')
    % TFN Workshop fit (1/m*atm)
    if (300 <= T)&&(T <= 2500) %Primer intervalo
        ap_g = 0.0;
    elseif (T > 2500) %Segundo intervalo
        ap_g = 0.0;
        warning('Calculo de ap fuera de rango')
    end
    
elseif strcmpi(tipo_gas,'O2')
    % TFN Workshop fit (1/m*atm)
    if (300 <= T)&&(T <= 2500) %Primer intervalo
        ap_g = 0.0;
    elseif (T > 2500) %Segundo intervalo
        ap_g = 0.0;
        warning('Calculo de ap fuera de rango')
    end
    
elseif strcmpi(tipo_gas,'CO2')
    % TFN Workshop fit (1/m*atm)
    if (300 <= T)&&(T <= 2500) %Primer intervalo
        c0 = 18.741; c1 = -121.310; c2 = 273.500; c3 = -194.050; c4 = 56.310; c5 = -5.8169;
        ap_g    = c0 + c1*(1000/T) + c2*(1000/T)^2 + c3*(1000/T)^3 + c4*(1000/T)^4 + c5*(1000/T)^5;
    elseif (T > 2500) %Segundo intervalo
        warning('Calculo de ap fuera de rango')
        c0 = 18.741; c1 = -121.310; c2 = 273.500; c3 = -194.050; c4 = 56.310; c5 = -5.8169;
        ap_g    = c0 + c1*(1000/2500) + c2*(1000/2500)^2 + c3*(1000/2500)^3 + c4*(1000/2500)^4 + c5*(1000/2500)^5;
    end
    
elseif strcmpi(tipo_gas,'H2O')
    % TFN Workshop fit (1/m*atm)
    if (300 <= T)&&(T <= 2500) %Primer intervalo
        c0      = -0.23093; c1 = -1.12390; c2 = 9.41530; c3 = -2.99880; c4 = 0.51382; c5 = -1.86840e-05;
        ap_g    = c0 + c1*(1000/T) + c2*(1000/T)^2 + c3*(1000/T)^3 + c4*(1000/T)^4 + c5*(1000/T)^5;
    elseif (T > 2500) %Segundo intervalo
        warning('Calculo de ap fuera de rango')
        c0      = -0.23093; c1 = -1.12390; c2 = 9.41530; c3 = -2.99880; c4 = 0.51382; c5 = -1.86840e-05;
        ap_g    = c0 + c1*(1000/2500) + c2*(1000/2500)^2 + c3*(1000/2500)^3 + c4*(1000/2500)^4 + c5*(1000/2500)^5;
    end
end

end

