function fes = FES_QX_Create(mesh, p)
    
    %Set type:
    fes.type        = 'QX';
    
    %Mesh:
    fes.mesh        = mesh;
    if mesh.nElems<2
%         error("There must be more than two elements")
    end
    
    %Discretization order:
    fes.p           = p;
    fes.nDof        = fes.mesh.nElems*(p+1);
    fes.ElemsDof    = reshape(1:fes.nDof, p+1, fes.mesh.nElems).';
    
    %Quadrature rule:
    fes.QuadRule    = GaussLegendre(2*p+1);
    
    %Functions to compute shape functions:
    fes.NCompute    = @(xi) PolyLegendre(xi, fes.p);
    fes.dNCompute   = @(xi) dPolyLegendre(xi, fes.p);
    
end

