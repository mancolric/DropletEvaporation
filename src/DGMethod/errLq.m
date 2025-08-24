function e_lq = errLq(u_h, NF, fes, u_theor, q)

    mesh        = fes.mesh;
    QuadRule    = fes.QuadRule;
    p           = fes.p;

    %Evaluate theoretical solution:
    x_qp        = PhysicalCoordinates(mesh, QuadRule.xi);
    u_qp        = u_theor(x_qp);

    %Numerical solution:
    uh_qp       = EvalSolution(u_h, fes, 1:mesh.nElems, QuadRule.xi);
    nVars       = length(uh_qp);

    %e_Lq^q |Omega| = int( 1/nVars * sum_I ((u_I^h - u_I^{theor})/NF_I)^q ) 
    omega_qp    = mesh.omega(x_qp);
    eq_qp       = zeros(mesh.nElems, QuadRule.N);
    for II=1:nVars
        eq_qp   = eq_qp + (abs(u_qp{II}-uh_qp{II})/NF(II)).^q;
    end
    eq_elems    = mesh.J .* ((eq_qp.*omega_qp)*QuadRule.w);

    Omega       = mesh.J.' * omega_qp * QuadRule.w;
    e_lq        = (sum(eq_elems)/nVars/Omega)^(1/q);

end