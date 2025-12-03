function value = devuelve_propiedad_gas(obj, propiedad, comb, estado, T )

for iD_combustible=1:size(obj.comp_inerts,2)
    if strcmpi(obj.comp_inerts{iD_combustible}, comb)
        indice_comb=iD_combustible;
        break
    end
end
MW  = obj.Inerts_mw(indice_comb);

if strcmpi(comb,'N2')

    % Propiedades N2(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')

        error("Liquid O2 ?")
        % Propiedades N2(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        % Entalpía de formación, 'Hf' [J/Kg]
        if strcmpi(propiedad,'Hf')
            value=0*1000/MW;
        else
            error("This property name does not exist")
        end

    end

elseif strcmpi(comb,'O2')

    % Propiedades O2(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')

        error("Liquid O2 ?")
        % Propiedades O2(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        % Entalpía de formación, 'Hf' [J/Kg]
        if strcmpi(propiedad,'Hf')
            value=0*1000/MW;
        else
            error("This property name does not exist")
        end

    end
elseif strcmpi(comb,'CO2')

    % Propiedades CO2(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')
        
        error("Liquid CO2 ?")
        % Propiedades CO2(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        % Entalpía de formación, 'Hf' [J/Kg]
        if strcmpi(propiedad,'Hf')
            value=-393.510*1000/MW;
        else
            error("This property name does not exist")
        end

    end
elseif strcmpi(comb,'H2O')

    % Propiedades H2O(líq.), todas en función de la Temperatura en K.
    if strcmpi(estado,'liquido')

        error("Liquid H2O ?")
        % Propiedades H2O(vap.), todas en función de la Temperatura en K.
    elseif strcmpi(estado,'vapor')
        % Entalpía de formación, 'Hf' [J/Kg] 
        if strcmpi(propiedad,'Hf')
            value=-241.826*1000/MW;
        else
            error("This property name does not exist")
        end

    end
end