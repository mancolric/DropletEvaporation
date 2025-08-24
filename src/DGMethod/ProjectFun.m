%Function to compute int(f * phi_i dOmega):
function b=ProjectFun(f, fes)

    mesh        = fes.mesh;
    p           = fes.p;
    QuadRule    = fes.QuadRule;
    
    %Location of quadrature nodes in the physical mesh:
    %Matrix is of size (Nx, NQp) (NQp is the number of quadrature nodes):
    x_qp        = PhysicalCoordinates(mesh, QuadRule.xi);
    omega_qp    = mesh.omega(x_qp);
    
    %Evaluate function at quadrature points:
    f_qp        = f(x_qp);
    nVars       = length(f_qp);
    
    %Load shape functions in reference element at quadrature nodes:
    %Matrix is of size (#quadrature nodes, p+1):
    phim        = fes.NCompute(QuadRule.xi);
    
    %For each element and degree of freedom, evaluate 
    %   sum_qp w_qp * f(x_gp) * phi_i * |J|:
    %The sum can be written as
    %   mesh.J .* ( f_qp * (w_qp.*phim[:,ii]) )
    %   [Nx,1] .* ( [Nx, Nqp] * ( [Nqp, 1] .* [Nqp, 1] ) )
    b           = zeros(nVars*fes.nDof, 1);
    for II=1:nVars
        for ii=1:p+1
            rows    = (II-1)*fes.nDof + fes.ElemsDof(:,ii);
            vals    = mesh.J .* ((omega_qp.*f_qp{II}) * (QuadRule.w.*phim(:,ii)));
            b(rows) = b(rows)+vals;
        end
    end
    
end