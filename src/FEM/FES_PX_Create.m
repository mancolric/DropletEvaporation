function fes = FES_PX_Create(mesh, p)
    
    %Set type:
    fes.type        = 'PX';
    
    %Mesh:
    fes.mesh        = mesh;
    if mesh.nElems<2
%         error("There must be more than two elements")
    end
    
    %Discretization order:
    fes.p           = p;
    fes.nDof        = fes.mesh.nElems*p+1;
    ElemsDof        = zeros(mesh.nElems, p+1);
    for iElem=1:mesh.nElems
        ElemsDof(iElem, :)  = (iElem-1)*p+1:(iElem-1)*p+p+1;
    end
    fes.ElemsDof    = ElemsDof;
    
    %Chebyshev nodes:
    fes.xi          = xi_Chebyshev(p);
    fes.NodesCoords = ChebyMesh(fes.mesh.x_faces, p);
    
    %Quadrature rule:
    fes.QuadRule    = GaussLegendre(2*p+1);
    
    %The Lagrange polynomials L_j(x) can be written as phi_k(x) a_kj, with
    %phi_k(x) the Legendre polynomials. The matrix a_kj can be computed
    %from
    %   delta_ij = L_j(x_i) = phi_k(x_i) a_kj,
    %where x_i are the nodes associated with the Lagrange polynomials.
    fes.Lag_Leg     = inv(PolyLegendre(fes.xi, p));
    
    %Functions to compute shape functions:
    fes.NCompute    = @(xi) PolyLag(xi, fes.p, fes.Lag_Leg);
    fes.dNCompute   = @(xi) dPolyLag(xi, fes.p, fes.Lag_Leg);
    
    %Matrix to assembly fast:
    fes.AsblyMat    = sparse(fes.ElemsDof(:), ...
                            (1:mesh.nElems*(p+1)).', ...
                            ones(mesh.nElems*(p+1),1), ...
                            fes.nDof, mesh.nElems*(p+1)); 
    
end

function Lm = PolyLag(xi, p, Lag_Leg)

    Lm      = PolyLegendre(xi, p)*Lag_Leg;
    
end

function dLm = dPolyLag(xi, p, Lag_Leg)

    dLm     = dPolyLegendre(xi, p)*Lag_Leg;
    
end
