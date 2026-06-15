classdef clase_IC
    %CLASE_IC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        diff_h
        nu
        cpFuel
        cpO2
        T_inf
        hfg
        kg
        rs
        Tb
    end
    
    methods
        function obj = clase_IC(gas, T_0, T_inf, x_lg)
            %CLASE_IC Construct an instance of this class
            %   Detailed explanation goes here
            % obj.nu      = 3.52;
            obj.nu      = gas.RStoi(1,2)* gas.species_mw(2) / (gas.RStoi(1,end) * gas.species_mw(end));     %Suponiendo un solo combustible
            obj.diff_h  = gas.gota.h_comb * gas.frac_masG{2};
            Cell_F          = cell(1, gas.nInerts + 1);
            Cell_F(1:end-1) = {0};  
            Cell_F{end}     = 1;
            obj.cpFuel  = calc_Cp(T_0, Cell_F, gas);
            % obj.cpO2    = calc_Cp(T_0, {0,1,0,0,0}, gas);
            obj.cpO2    = 1175;
            obj.T_inf   = T_inf;
            obj.hfg     = gas.Lv;
            % obj.kg      = calc_D_T(T_0, {0,0,0,0,1}, gas);
            obj.kg      = MixtureRules('k_gas', T_0, zeros(gas.nInerts, 1), 1, gas.matrix, gas.gota, gas.comp_inerts, gas.P);
            obj.rs      = x_lg;
            obj.Tb      = gas.gota.Tb;
        end
        
    end
end

