function [J,Jnum] = FluxSourceRestriction_CheckJacobian(model, t, u, fes, NF)

    %Evaluate r.h.s. and Jacobian:
    [~,J]       = FluxSourceRestriction(model, t, u, fes, true);
    
    %Numerical Jacobian:
    delta       = 1e-8;
    Jnum        = zeros(size(J));
    uv          = CellToVector(u);
    nVars       = length(u);
    for iVar=1:length(u)
        for iDof=1:fes.nDof
            ii          = (iVar-1)*fes.nDof + iDof;
            upert       = uv;
            upert(ii)   = upert(ii)-NF(iVar)*delta;
            [b1,~]      = FluxSourceRestriction(model, t, VectorToCell(upert,nVars), fes, true);
            upert       = uv;
            upert(ii)   = upert(ii)+NF(iVar)*delta;
            [b2,~]      = FluxSourceRestriction(model, t, VectorToCell(upert,nVars), fes, true);
            Jnum(:,ii)  = (b2-b1)/(2*NF(iVar)*delta);
        end
    end
    
    %Print info:
    disp(['||J-J_num||=',num2str(norm(Jnum-J))])
    disp(['||J_num||=',num2str(norm(Jnum))])
    
end
    