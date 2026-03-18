%Function to compute int(omega*phi_i*phi_j dOmega):
function M = MassMatrix(fes)

    mesh        = fes.mesh;
    p           = fes.p;
    QuadRule    = fes.QuadRule;
    
    %Location of quadrature nodes in the physical mesh:
    %Matrix is of size (Nx, NQp) (NQp is the number of quadrature nodes):
    x_qp        = PhysicalCoordinates(mesh, QuadRule.xi);
    
    %Evaluate function omega at quadrature points:
    omega_qp    = mesh.omega(x_qp);
    
    %Load shape functions in reference element at quadrature nodes:
    %Matrix is of size (#quadrature nodes, p+1):
    phim        = fes.NCompute(QuadRule.xi);
    
    %For each element and degree of freedom, evaluate 
    %   sum_qp w_qp * omega(x_gp) * phi_i * phi_j |J|.
    rows        = zeros(mesh.nElems, p+1, p+1);
    cols        = zeros(mesh.nElems, p+1, p+1);
    vals        = zeros(mesh.nElems, p+1, p+1);
    for ii=1:p+1
        for jj=1:p+1
            
            %For element K, local node "i" corresponds to global node
            %ElemsDof(K,ii)
            rows(:,ii,jj)   = fes.ElemsDof(:,ii);
            cols(:,ii,jj)   = fes.ElemsDof(:,jj);
            
            %sum_qp w_qp * omega(x_qp) * phi_i * phi_j |J| can be written as
            %   mesh.J .* ( omega_qp * (w_qp.*phim[:,ii].*phim[:,jj]) )
            %   [Nx,1] .* ( [Nx, Nqp] * ( [Nqp, 1] .* [Nqp, 1] .* [Nqp, 1] )
            vals(:,ii,jj)   = mesh.J .* (omega_qp * (QuadRule.w .* (phim(:,ii).*phim(:,jj))) );
        
        end
    end
    
    M       = sparse(rows(:), cols(:), vals(:));
    
end