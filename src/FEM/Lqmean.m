%Compute Lq mean of a solution:
function umean = Lqmean(u, fes, q)

    mesh        = fes.mesh;
    QuadRule    = fes.QuadRule;
    
    x_qp        = PhysicalCoordinates(mesh, QuadRule.xi);
    omega_qp    = mesh.omega(x_qp);
    
    %Numerical solution:
    u_qp        = EvalSolution(u, fes, 1:mesh.nElems, QuadRule.xi);
    nVars       = length(u_qp);
    umean       = zeros(nVars, 1);
    
    %Return maximum values if q=Inf:
    if q==Inf
        for II=1:nVars
            umean(II)   = norm(u_qp{II},Inf);
        end
        return
    end
    
    %Compute domain volume:
    Omega       = mesh.J.' * omega_qp * QuadRule.w;
    
    %Loop variables and compute umean = (int(|u|^q dOmega)/Omega)^(1/q):
    for II=1:nVars
        umean(II)   = (mesh.J.' * (abs(u_qp{II}).^q .* omega_qp) * QuadRule.w / Omega)^(1/q);
    end
    
end
