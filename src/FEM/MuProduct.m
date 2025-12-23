function [Mu, dMu_du, dMu_dx] = MuProduct(fes, u, ComputeJ)

    %Extract values:
    mesh        = fes.mesh;
    p           = fes.p;
    QuadRule    = fes.QuadRule;
    nVars       = length(u);
    
    %Load shape functions in reference element at quadrature nodes:
    %Matrix is of size (#quadrature nodes, p+1):
    phim        = fes.NCompute(QuadRule.xi);
    phi1m       = cat(2, 0.5*(1.0-QuadRule.xi(:)), 0.5*(1.0+QuadRule.xi(:)));
    
    %For each element [a,b] and degree of freedom, evaluate 
    %   sum_qp w_qp * 4*pi* r^2 * u * phi_i * |J| = 
    %   sum_qp |J| (r^2 u_qp) (w_qp phi_i)
    %           [nElems,1] [nElems, nqp] [nqp,1] 
    Mu          = zeros(nVars*fes.nDof,1);
    r2_qp       = 4*pi*(mesh.x_faces(1:end-1)*phi1m(:,1).' + ...
                    mesh.x_faces(2:end)*phi1m(:,2).').^2;   %[nElems, nqp]
    for II=1:nVars
        Mu_ElemsDof = zeros(mesh.nElems, p+1);
        u_qp        = u{II}(fes.ElemsDof) * phim.';
        for ii=1:p+1
            Mu_ElemsDof(:,ii)   = mesh.J .* ( (r2_qp.*u_qp) * (QuadRule.w .* phim(:,ii)) );
        end
        dof         = (II-1)*fes.nDof+1:II*fes.nDof;
        Mu(dof)     = fes.AsblyMat*Mu_ElemsDof(:);
    end
    
    %Evaluate derivatives:
    if ComputeJ
        
        %Derivatives w.r.t. to u_j:
        %   sum_qp w_qp * 4*pi* r^2 * phi_j * phi_i * |J|
        %   sum_qp |J| (4*pi* r^2) * (w_qp * phi_j * phi_i)
        dMu_du.iv   = zeros(mesh.nElems, p+1, p+1, nVars);
        dMu_du.jv   = zeros(mesh.nElems, p+1, p+1, nVars);
        dMu_du.sv   = zeros(mesh.nElems, p+1, p+1, nVars);
        for II=1:nVars
            cum_dof = (II-1)*fes.nDof;
            for jj=1:p+1
                for ii=1:p+1
                    dMu_du.iv(:,ii,jj,II)   = cum_dof + fes.ElemsDof(:,ii);
                    dMu_du.jv(:,ii,jj,II)   = cum_dof + fes.ElemsDof(:,jj);
                    dMu_du.sv(:,ii,jj,II)   = mesh.J .* ( r2_qp * (QuadRule.w .* phim(:,ii) .* phim(:,jj) ) );
                end
            end
        end
        
        %Derivatives w.r.t. to x_j:
        %   sum_qp w_qp * 4*pi* 2r L_j * u * phi_i * |J| + 
        %   sum_qp w_qp * 4*pi r^2 u phi_i * dL_j/dxi
        %
        %   sum_qp |J| (8*pi* r * u) * (w_qp * L_j * phi_i) + 
        %   sum_qp (4*pi r^2 u) * (w_qp * dL_j/dxi * phi_i)
        r_qp        = 4*pi*(mesh.x_faces(1:end-1)*phi1m(:,1).' + ...
                        mesh.x_faces(2:end)*phi1m(:,2).');   %[nElems, nqp]
        dMu_dx.iv   = zeros(mesh.nElems, p+1, 2, nVars);
        dMu_dx.jv   = zeros(mesh.nElems, p+1, 2, nVars);
        dMu_dx.sv   = zeros(mesh.nElems, p+1, 2, nVars);
        ElemsNodes  = cat(2, (1:mesh.nElems).', (2:mesh.nElems+1).');
        dphi1_dxi   = [ -0.5, +0.5 ];
        for II=1:nVars
            cum_dof     = (II-1)*fes.nDof;
            u_qp        = u{II}(fes.ElemsDof) * phim.';
            for jj=1:2
                for ii=1:p+1
                    dMu_dx.iv(:,ii,jj,II)   = cum_dof + fes.ElemsDof(:,ii);
                    dMu_dx.jv(:,ii,jj,II)   = ElemsNodes(:,jj);
                    dMu_dx.sv(:,ii,jj,II)   = ...
                        mesh.J .* ( (2.0 .* r_qp .* u_qp) * (QuadRule.w .* phim(:,ii) .* phi1m(:,jj) ) ) + ...
                        (r2_qp .* u_qp) * (QuadRule.w .* phim(:,ii) .* dphi1_dxi(jj));
                end
            end
        end
        
        %Reshape:
        dMu_du.iv   = dMu_du.iv(:);
        dMu_du.jv   = dMu_du.jv(:);
        dMu_du.sv   = dMu_du.sv(:);
        dMu_dx.iv   = dMu_dx.iv(:);
        dMu_dx.jv   = dMu_dx.jv(:);
        dMu_dx.sv   = dMu_dx.sv(:);
        
    else
        
        dMu_du      = NaN;
        dMu_dx      = NaN;
        
    end     
    
end