function [Mm, usol] = InitialCondition(fes, u0fun)

    %Compute mass matrix for one variable:
    Mm          = MassMatrix(fes);
    
    %Project all variables:
    b           = ProjectFun(u0fun, fes);
    
    %Solve each variable:
    nVars       = length(b)/fes.nDof;
    usol        = cell(nVars,1);
    for II=1:nVars
        usol{II}    = Mm\b((II-1)*fes.nDof+1:II*fes.nDof);
    end
    
end