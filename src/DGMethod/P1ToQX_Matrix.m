function P1ToQXmat = P1ToQX_Matrix(fes)

    %If a0, a1 are the Legendre coefficients, and u1 and u2 are the values
    %at the Lagrangian nodes, 
    %   a0      = 0.5*(u1+u2)
    %   a1      = 0.5*(u2-u1);
    
    nElems      = fes.mesh.nElems;
    p           = fes.p;
    
    %Influence of u1 on a0:
    i01v        = 1:(p+1):((nElems-1)*(p+1)+1);
    j01v        = 1:nElems;
    s01v        = 0.5 + 0.0*i01v;
    %Influence of u1 on a1:
    i11v        = 2:(p+1):((nElems-1)*(p+1)+2);
    j11v        = 1:nElems;
    s11v        = - 0.5 + 0.0*i11v;
    %Influence of u2 on a0:
    i02v        = 1:(p+1):((nElems-1)*(p+1)+1);
    j02v        = 2:nElems+1;
    s02v        = 0.5 + 0.0*i01v;
    %Influence of u2 on a1:
    i12v        = 2:(p+1):((nElems-1)*(p+1)+2);
    j12v        = 2:nElems+1;
    s12v        = 0.5 + 0.0*i11v;
    
    %Construct sparse matrix:
    P1ToQXmat   = sparse(   cat(2, i01v, i11v, i02v, i12v), ...
                            cat(2, j01v, j11v, j02v, j12v), ...
                            cat(2, s01v, s11v, s02v, s12v), ...
                            nElems*(p+1), nElems+1) ;
                        
end
 