function unodes = P1ToPX(uverts, fes)
    
    nElems      = fes.mesh.nElems;
    p           = fes.p;
    unodes      = zeros(nElems*p+1,1);
    
    %Lagrangian shape functions of P1 at Chebyshev nodes of PX:
    N1v         = 0.5*(1.0-fes.xi);
    N2v         = 0.5*(1.0+fes.xi);
    u1v         = uverts(1:end-1);
    u2v         = uverts(2:end);
    %Loop inner nodes:
    for iDof=2:p
        dof                 = (0:nElems-1)*p + iDof;
        unodes(dof)         = u1v*N1v(iDof) + u2v*N2v(iDof);
    end
    %Vertices:
    unodes(1:p:nElems*p+1)  = uverts;
    
end
 