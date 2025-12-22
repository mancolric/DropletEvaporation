function u_Lq = LqNorm(u_h, NF, fes, q)

    mesh        = fes.mesh;
    QuadRule    = fes.QuadRule;

    %Numerical solution:
    x_qp        = PhysicalCoordinates(mesh, QuadRule.xi);
    uh_qp       = EvalSolution(u_h, fes, 1:mesh.nElems, QuadRule.xi);
    nVars       = length(uh_qp);

    %Return maximum value if q=Inf
    if q==Inf
        u_Lq    = 0.0;
        for II=1:nVars
            u_Lq    = max(u_Lq, norm(uh_qp{II}(:), Inf)/NF(II));
        end
        return
    end
    
    %nVars u_Lq^q |Omega| = int( sum_I (u_I^h/NF_I)^q ) 
    omega_qp    = mesh.omega(x_qp);
    uq_qp       = zeros(mesh.nElems, QuadRule.N);
    for II=1:nVars
        uq_qp   = uq_qp + (abs(uh_qp{II})/NF(II)).^q;
    end
    uq_elems    = mesh.J .* ((uq_qp.*omega_qp)*QuadRule.w);

    Omega       = mesh.J.' * omega_qp * QuadRule.w;
    u_Lq        = (sum(uq_elems)/nVars/Omega)^(1/q);

end
