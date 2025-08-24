function unodes = QXToP1_Matrix(uLeg, fes)

    error('Unfinished')
    
    uLeg                = uLeg(fes.ElemsDof);
    uL                  = uLeg(:,1)-uLeg(:,2);
    uR                  = uLeg(:,1)+uLeg(:,2);
    unodes              = zeros(fes.mesh.nElems+1,1);
    unodes(1)           = uL(1);
    unodes(2:end-1)     = 0.5*(uL(2:end)+uR(1:end-1));
    unodes(end)         = uR(end);
    
end
 