function [Unknowns] = InitialCond(obj, ctes_ini)
%Unknowns = zeros(5,1); %[mdot, Ts, Tf, Yfs, rf]


function [l, varargout] = Residue(ctes, bool)
    l           = zeros(length(ctes),1);
    B_oq        =  ((obj.diff_h / obj.nu)+(obj.cpg*(obj.T_inf - ctes(2))))/obj.hfg;
    if nargout==1
        l(1)            = ctes(1) - 4*pi*obj.kg*obj.rs*log(1+B_oq)/obj.cpg;
        l(2)            = obj.hfg*(obj.nu*B_oq-1)/(obj.cpg*(1+obj.nu)) + ctes(2) - ctes(3);
        l(3)            = ctes(5) - obj.rs * log(1+B_oq)/log((obj.nu + 1)/obj.nu);
        l(4)            = ctes(4) - (B_oq - (1/obj.nu))/(B_oq + 1);
        l(5)            = obj.Tb - 20 - ctes(2);
    else
        l(1)            = ctes(1) - 4*pi*obj.kg*obj.rs*log(1+B_oq)/obj.cpg;
        l(2)            = obj.hfg*(obj.nu*B_oq-1)/(obj.cpg*(1+obj.nu)) + ctes(2) - ctes(3);
        l(3)            = ctes(5) - obj.rs * log(1+B_oq)/log((obj.nu + 1)/obj.nu);
        l(4)            = ctes(4) - (B_oq - (1/obj.nu))/(B_oq + 1);
        l(5)            = obj.Tb - 20 - ctes(2);

        varargout{1}    = ComputeJacobian(@(ctes2) Residue(ctes2), ctes, 1e-6);
    end
end

[Unknowns, Iteration, flag]    = NewtonRaphsonLS(@Residue, ctes_ini, 0.0, 1e-8, 100, 100);

if flag<0
    error("En el cálculo de perfil inicial líquido")
end

end

