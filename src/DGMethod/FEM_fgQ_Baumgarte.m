function [F, dF_dU, dF_dq1, dF_dqN, Deltat_CFL] = FEM_fgQ_Baumgarte(model, t, usol, fes, ComputeJ, ...
    Mm, NF, tau)

    %NOTE: This function is only valid for DG.
    
    mesh        = fes.mesh;
    p           = fes.p;
    QuadRule    = fes.QuadRule;
    nElems      = mesh.nElems;
    
    %Allocate rhs vector:
    nVars       = length(usol);
    F           = zeros((model.nDiff+model.nAlg)*fes.nDof, 1);
    
    %----------------------------------------------------------------------
    %DOMAIN INTEGRALS:
    
    %Location of quadrature nodes in the physical mesh:
    %Matrix is of size (nElems, NQp) (NQp is the number of quadrature nodes):
    x_qp        = PhysicalCoordinates(mesh, QuadRule.xi);
    omega_qp    = mesh.omega(x_qp);
    Jinv_qp     = (1.0./mesh.J) * ones(1,QuadRule.N);
    womegaJ_qp  = omega_qp .* (mesh.J*QuadRule.w.');
    
    %Evaluate solution at quadrature points:
    u_qp        = EvalSolution(usol, fes, 1:mesh.nElems, QuadRule.xi); %cell(nVars,1)
    du_dx_qp    = EvalSolution_dx(usol, fes, 1:mesh.nElems, QuadRule.xi); %cell(nVars,1)
    
    %Evaluate flux and source terms at quadrature points:
    [f_qp, df_du_qp, df_dgradu_qp, ...
        Q_qp, dQ_du_qp, dQ_dgradu_qp, ...
        g_qp, dg_du_qp, dg_dgradu_qp, lambdav] ...
                = model.fQg(model, t, x_qp, u_qp, du_dx_qp, ComputeJ);
    
    %Load shape functions in reference element at quadrature nodes:
    %Matrix is of size (#quadrature nodes, p+1):
    psim        = fes.NCompute(QuadRule.xi);
    dpsim       = fes.dNCompute(QuadRule.xi);
    
    %Contribution of flux and source terms in the domain:
    %   int(Q_I*psi_alpha) + int(f_I*dpsi_alpha/dx)
    for II=1:model.nDiff
        F((II-1)*fes.nDof+fes.ElemsDof)     = (f_qp{II}.*Jinv_qp.*womegaJ_qp) * dpsim + ...
                                                (Q_qp{II}.*womegaJ_qp) * psim;
    end
    
    %Compute db_(I,alpha)/du_(J,beta). 
    if ComputeJ
        
        %Allocate memory for jacobian elements. For each (II,JJ) block, we have
        %nElems*(p+1)^2 nonzeros:
        iv      = zeros(mesh.nElems, (p+1)*(p+1), model.nDiff, nVars);
        jv      = zeros(size(iv));
        sv      = zeros(size(iv));
        
        %Products of shape functions:
        UpsilonNN       = UpsilonMat(psim, psim);
        UpsilonGN       = UpsilonMat(dpsim, psim);
        UpsilonNG       = UpsilonMat(psim, dpsim);
        UpsilonGG       = UpsilonMat(dpsim, dpsim);
       
        %Loop blocks:
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof, fes.ElemsDof);
        for II=1:model.nDiff
            for JJ=1:nVars
                
                %Save rows and columns of nnz:
                iv(:,:,II,JJ)   = imat + (II-1)*fes.nDof;
                jv(:,:,II,JJ)   = jmat + (JJ-1)*fes.nDof;
                
                %Four contributions:
                %   int df_I/du_J dxi/dx dpsi_alpha/dxi psi_beta
                %   int df_I/d(du_J) (dxi/dx)^2 dpsi_alpha/dxi dpsi_beta/dxi
                %   int dQ_I/du_J psi_alpha psi_beta
                %   int dQ_I/d(du_J) dxi/dx psi_alpha dpsi_beta
                sv(:,:,II,JJ)   = (df_du_qp{II,JJ} .* Jinv_qp .* womegaJ_qp) * UpsilonGN + ...
                                    (df_dgradu_qp{II,JJ} .* Jinv_qp.^2 .* womegaJ_qp) * UpsilonGG + ...
                                    (dQ_du_qp{II,JJ} .* womegaJ_qp) * UpsilonNN + ...
                                    (dQ_dgradu_qp{II,JJ} .* Jinv_qp .* womegaJ_qp) * UpsilonNG;
            end
        end
        
    end
    
    %----------------------------------------------------------------------
    %FLUXES AT INTERNAL FACES:
    
    %Solution at west and east faces:
    u_w         = EvalSolution(usol, fes, 2:nElems, [-1.0]); %cell(nVars,1)
    du_dx_w     = EvalSolution_dx(usol, fes, 2:nElems, [-1.0]); %cell(nVars,1)
    u_e         = EvalSolution(usol, fes, 1:nElems-1, [1.0]); %cell(nVars,1)
    du_dx_e     = EvalSolution_dx(usol, fes, 1:nElems-1, [1.0]); %cell(nVars,1)
    psim_w      = fes.NCompute([-1.0]);
    psim_e      = fes.NCompute([1.0]);
    dpsim_w     = fes.dNCompute([-1.0]);
    dpsim_e     = fes.dNCompute([1.0]);
    Jinv_elems  = 1.0./mesh.J;
    omega_faces = mesh.omega(mesh.x_faces(2:end-1));
    hp_elems    = diff(mesh.x_faces)/fes.p;
    hp_faces    = 2.0./(1.0./hp_elems(1:nElems-1)+1.0./hp_elems(2:nElems));

    %Fluxes at the nElems-1 internal faces:
    [fhat, dfhat_du_e, dfhat_dgradu_e, ...
        dfhat_du_w, dfhat_dgradu_w]     = model.ftilde(model, t, mesh.x_faces(2:end-1), ...
                                            u_e, du_dx_e, u_w, du_dx_w, hp_faces, ...
                                            ComputeJ);
                                        
    %Fluxes are leaving the elements 1:nElems-1 at the east faces:
    for II=1:model.nDiff
        dof         = (II-1)*fes.nDof+fes.ElemsDof(1:nElems-1,:);
        F(dof)      = F(dof) - (fhat{II}.*omega_faces) * psim_e;
    end
    if ComputeJ
        
        %Variations of leaving fluxes w.r.t. variations at left:
        iv_ee           = zeros(nElems-1, (p+1)*(p+1), model.nDiff, nVars);
        jv_ee           = zeros(size(iv_ee));
        sv_ee           = zeros(size(iv_ee));
        UpsilonNN       = UpsilonMat(psim_e, psim_e);
        UpsilonNG       = UpsilonMat(psim_e, dpsim_e);
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof(1:nElems-1,:), fes.ElemsDof(1:nElems-1,:));
        for II=1:model.nDiff
            for JJ=1:nVars
                iv_ee(:,:,II,JJ)    = imat + (II-1)*fes.nDof;
                jv_ee(:,:,II,JJ)    = jmat + (JJ-1)*fes.nDof;
                sv_ee(:,:,II,JJ)    = - (dfhat_du_e{II,JJ}.*omega_faces) * UpsilonNN - ...
                                        (dfhat_dgradu_e{II,JJ}.*Jinv_elems(1:nElems-1).*omega_faces) * UpsilonNG;
            end
        end
       
        %Variations of leaving fluxes w.r.t. variations at the right:
        iv_ew           = zeros(nElems-1, (p+1)*(p+1), model.nDiff, nVars);
        jv_ew           = zeros(size(iv_ew));
        sv_ew           = zeros(size(iv_ew));
        UpsilonNN       = UpsilonMat(psim_e, psim_w);
        UpsilonNG       = UpsilonMat(psim_e, dpsim_w);
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof(1:nElems-1,:), fes.ElemsDof(2:nElems,:));
        for II=1:model.nDiff
            for JJ=1:nVars
                iv_ew(:,:,II,JJ)    = imat + (II-1)*fes.nDof;
                jv_ew(:,:,II,JJ)    = jmat + (JJ-1)*fes.nDof;
                sv_ew(:,:,II,JJ)    = - (dfhat_du_w{II,JJ}.*omega_faces) * UpsilonNN - ...
                                        (dfhat_dgradu_w{II,JJ}.*Jinv_elems(2:nElems).*omega_faces) * UpsilonNG;
            end
        end
        
    end
    
    %Fluxes are entering the elements 2:nElems at the west faces:
    for II=1:model.nDiff
        dof         = (II-1)*fes.nDof+fes.ElemsDof(2:nElems,:);
        F(dof)      = F(dof) + (fhat{II}.*omega_faces) * psim_w;
    end
    if ComputeJ
        
        %Variations of entering fluxes w.r.t. variations at left:
        iv_we           = zeros(nElems-1, (p+1)*(p+1), model.nDiff, nVars);
        jv_we           = zeros(size(iv_we));
        sv_we           = zeros(size(iv_we));
        UpsilonNN       = UpsilonMat(psim_w, psim_e);
        UpsilonNG       = UpsilonMat(psim_w, dpsim_e);
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof(2:nElems,:), fes.ElemsDof(1:nElems-1,:));
        for II=1:model.nDiff
            for JJ=1:nVars
                iv_we(:,:,II,JJ)    = imat + (II-1)*fes.nDof;
                jv_we(:,:,II,JJ)    = jmat + (JJ-1)*fes.nDof;
                sv_we(:,:,II,JJ)    = + (dfhat_du_e{II,JJ}.*omega_faces) * UpsilonNN + ...
                                        (dfhat_dgradu_e{II,JJ}.*Jinv_elems(1:nElems-1).*omega_faces) * UpsilonNG;
            end
        end
       
        %Variations of entering fluxes w.r.t. variations at the right:
        iv_ww           = zeros(nElems-1, (p+1)*(p+1), model.nDiff, nVars);
        jv_ww           = zeros(size(iv_ww));
        sv_ww           = zeros(size(iv_ww));
        UpsilonNN       = UpsilonMat(psim_w, psim_w);
        UpsilonNG       = UpsilonMat(psim_w, dpsim_w);
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof(2:nElems,:), fes.ElemsDof(2:nElems,:));
        for II=1:model.nDiff
            for JJ=1:nVars
                iv_ww(:,:,II,JJ)    = imat + (II-1)*fes.nDof;
                jv_ww(:,:,II,JJ)    = jmat + (JJ-1)*fes.nDof;
                sv_ww(:,:,II,JJ)    = + (dfhat_du_w{II,JJ}.*omega_faces) * UpsilonNN + ...
                                        (dfhat_dgradu_w{II,JJ}.*Jinv_elems(2:nElems).*omega_faces) * UpsilonNG;
            end
        end
        
    end
    
    %----------------------------------------------------------------------
    %FLUXES AT BOUNDARY FACES (I.E., IMPOSITION OF B.C.):
    
    %Left boundary:
    u1              = EvalSolution(usol, fes, [1], [-1.0]); %cell(nVars,1)
    du1_dx          = EvalSolution_dx(usol, fes, [1], [-1.0]); %cell(nVars,1)
    omega_face1     = mesh.omega(mesh.x_faces(1));
    
    %Flux at left boundary is entering the first element:
    [f1, df1_du1, df1_du1_dx, df1_dq1] = ...
        model.ftilde1(model, t, mesh.x_faces(1), ...
                        u1, du1_dx, hp_elems(1), ComputeJ);
    for II=1:model.nDiff
        dof     = fes.ElemsDof(1,:) + (II-1)*fes.nDof;
        F(dof)  = F(dof) + reshape((f1{II}.*omega_face1) * psim_w, [], 1); %reshape is necessary because it is only one element
    end
    if ComputeJ
        %Variations of incoming flux at face 1 w.r.t. variations in
        %solution at face 1:
        iv_11           = zeros(1, (p+1)*(p+1), model.nDiff, nVars);
        jv_11           = zeros(size(iv_11));
        sv_11           = zeros(size(iv_11));
        UpsilonNN       = UpsilonMat(psim_w, psim_w);
        UpsilonNG       = UpsilonMat(psim_w, dpsim_w);
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof(1,:), fes.ElemsDof(1,:));
        for II=1:model.nDiff
            for JJ=1:nVars
                iv_11(:,:,II,JJ)    = imat + (II-1)*fes.nDof;
                jv_11(:,:,II,JJ)    = jmat + (JJ-1)*fes.nDof;
                sv_11(:,:,II,JJ)    = (df1_du1{II,JJ} .* omega_face1) * ...
                                        UpsilonNN + ...
                                        (df1_du1_dx{II,JJ} .* Jinv_elems(1) .* omega_face1) * ...
                                        UpsilonNG;
            end
        end
        
        %Variations of flux at face 1 w.r.t. variations in the parameters
        %of the boundary condition:
        nParam          = size(df1_dq1,2);
        dF_dq1.iv       = zeros(p+1, model.nDiff, nParam);
        dF_dq1.jv       = zeros(p+1, model.nDiff, nParam);
        dF_dq1.sv       = zeros(p+1, model.nDiff, nParam);
        for II=1:model.nDiff
            for JJ=1:nParam
                dF_dq1.iv(:,II,JJ)    = fes.ElemsDof(1,:) + (II-1)*fes.nDof;
                dF_dq1.jv(:,II,JJ)    = JJ;
                dF_dq1.sv(:,II,JJ)    = reshape((df1_dq1{II,JJ}.*omega_face1) * psim_w, [], 1);
            end
        end
        dF_dq1.iv       = dF_dq1.iv(:);
        dF_dq1.jv       = dF_dq1.jv(:);
        dF_dq1.sv       = dF_dq1.sv(:);
        
    else
        
        dF_dq1          = NaN;
                
    end
    
    %Right boundary:
    uN              = EvalSolution(usol, fes, [nElems], [1.0]); %cell(nVars,1)
    duN_dx          = EvalSolution_dx(usol, fes, [nElems], [1.0]); %cell(nVars,1)
    omega_faceN     = mesh.omega(mesh.x_faces(nElems+1));
    
    %Flux at right boundary is leaving the last element:
    [fN, dfN_duN, dfN_duN_dx, dfN_dqN] = ...
        model.ftildeN(model, t, mesh.x_faces(end), ...
                        uN, duN_dx, hp_elems(nElems), ComputeJ);
    for II=1:model.nDiff
        dof     = fes.ElemsDof(nElems,:) + (II-1)*fes.nDof;
        F(dof)  = F(dof) - reshape((fN{II}.*omega_faceN) * psim_e, [], 1); %reshape is necessary because it is only one element
    end
    if ComputeJ
        
        %Variations of flux at face N w.r.t. variations in
        %solution at face N:
        iv_22           = zeros(1, (p+1)*(p+1), model.nDiff, nVars);
        jv_22           = zeros(size(iv_22));
        sv_22           = zeros(size(iv_22));
        UpsilonNN       = UpsilonMat(psim_e, psim_e);
        UpsilonNG       = UpsilonMat(psim_e, dpsim_e);
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof(nElems,:), fes.ElemsDof(nElems,:));
        for II=1:model.nDiff
            for JJ=1:nVars
                iv_22(:,:,II,JJ)    = imat + (II-1)*fes.nDof;
                jv_22(:,:,II,JJ)    = jmat + (JJ-1)*fes.nDof;
                sv_22(:,:,II,JJ)    = - (dfN_duN{II,JJ} .* omega_faceN) * ...
                                        UpsilonNN - ...
                                        (dfN_duN_dx{II,JJ} .* Jinv_elems(nElems) .* omega_faceN) * ...
                                        UpsilonNG;
            end
        end
        
        %Variations of flux at face N w.r.t. variations in the parameters
        %of the boundary condition:
        nParam          = size(dfN_dqN,2);
        dF_dqN.iv       = zeros(p+1, model.nDiff, nParam);
        dF_dqN.jv       = zeros(p+1, model.nDiff, nParam);
        dF_dqN.sv       = zeros(p+1, model.nDiff, nParam);
        for II=1:model.nDiff
            for JJ=1:nParam
                dF_dqN.iv(:,II,JJ)    = fes.ElemsDof(nElems,:) + (II-1)*fes.nDof;
                dF_dqN.jv(:,II,JJ)    = JJ;
                dF_dqN.sv(:,II,JJ)    = -reshape((dfN_dqN{II,JJ}.*omega_faceN) * psim_e, [], 1);
            end
        end
        dF_dqN.iv       = dF_dqN.iv(:);
        dF_dqN.jv       = dF_dqN.jv(:);
        dF_dqN.sv       = dF_dqN.sv(:);
        
    else
        
        dF_dqN          = NaN;
        
    end
    
    %----------------------------------------------------------------------
    %RESTRICTION:
    
    %Instead of considering the restriction G(U)=0, which yields a DAE2, 
    %we consider
    %   H(U,V)  = tau * dG/dU(U)*Udot(U,V) + G(U) = 
    %           = tau * dG/dU(U)*M^{-1}*F(U,V) + G(U) = 0
    %The term dG/dU(U)*M^{-1}*F(U,V) is estimated via high-order finite
    %differences.
    
    %Compute G(U) and derivatives:
    G           = zeros(model.nAlg*fes.nDof, 1);
    for II=1:model.nAlg
        dof     = (II-1)*fes.nDof + fes.ElemsDof;
        G(dof)  = (g_qp{II}.*womegaJ_qp) * psim;
    end
    if ComputeJ
        
        %Allocate memory:
        iv_g    = zeros(mesh.nElems, (p+1)*(p+1), model.nAlg, model.nDiff);
        jv_g    = zeros(size(iv_g));
        sv_g    = zeros(size(iv_g));
        
        %Products of shape functions:
        UpsilonNN       = UpsilonMat(psim, psim);
        UpsilonNG       = UpsilonMat(psim, dpsim);
        
        %int(dg/du psi psi + dg/d(du/dx) dpsi/dx psi):
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof, fes.ElemsDof);
        for II=1:model.nAlg
            for JJ=1:model.nDiff
                iv_g(:,:,II,JJ)     = imat + (II-1)*fes.nDof;
                jv_g(:,:,II,JJ)     = jmat + (JJ-1)*fes.nDof;
                sv_g(:,:,II,JJ)     = (dg_du_qp{II,JJ}.*womegaJ_qp) * UpsilonNN + ...
                                        (dg_dgradu_qp{II,JJ}.*Jinv_qp.*womegaJ_qp) * UpsilonNG;
            end
        end
        
    end
    
    %Compute Udot = M^{-1} F(U,V):
    Udot                = cell(model.nDiff, 1);
    Mm_F                = LUFactorization(Mm);
    Udot_rnorm          = zeros(model.nDiff, 1);
    for II=1:model.nDiff
        dof             = (II-1)*fes.nDof + (1:fes.nDof);
        Udot{II}        = LUSolve(Mm_F, F(dof));
        Udot_rnorm(II)  = norm(Udot{II},Inf)/NF(II);
    end
    udot_qp         = EvalSolution(Udot, fes, 1:mesh.nElems, QuadRule.xi);
    dudot_dx_qp     = EvalSolution_dx(Udot, fes, 1:mesh.nElems, QuadRule.xi);
    
    %Scale Udot's in such a way that the maximum relative norm is one:
    Udot_scal       = max(Udot_rnorm);
    
    %Derivation nodes and weights:
    epsilon         = 1e-4;
    xi_deriv        = [-2.0, -1.0, 0.0, 1.0, 2.0]*epsilon;
    w_deriv         = [1, -8, 0, 8, -1]/(12*epsilon);
    
    %Compute Gdot = dG/dU*Udot = ||Udot|| sum w_i G(U+xi_i*Udot/||Udot||):
    Gdot            = zeros(model.nAlg*fes.nDof, 1);
    Gpert           = zeros(model.nAlg*fes.nDof, 1);
    for ideriv=1:length(w_deriv)
        
        %Compute G(U+xi(ideriv)*Udot):
        if xi_deriv(ideriv)~=0.0
            
            %Compute U+xi_i*Udot:
            upert_qp        = cell(model.nVars, 1);
            dupert_dx_qp    = cell(model.nVars, 1);
            for II=1:model.nDiff
                upert_qp{II}        = u_qp{II} + xi_deriv(ideriv)*(udot_qp{II}/Udot_scal);
                dupert_dx_qp{II}    = du_dx_qp{II} + xi_deriv(ideriv)*(dudot_dx_qp{II}/Udot_scal);
            end
            for II=model.nDiff+1:model.nVars
                upert_qp{II}        = u_qp{II}*NaN;
                dupert_dx_qp{II}    = du_dx_qp{II}*NaN;
            end

            %Evaluate restriction g(upert) and corresponding term G(Upert)
%             g_qp            = model.g(model, t, x_qp, upert_qp, du_dx_pert, false);
            [~, ~, ~, ~, ~, ~, g_qp, ~, ~, ~] = ...
                model.fQg(model, t, x_qp, upert_qp, dupert_dx_qp, false);
            for II=1:model.nAlg
                dof         = (II-1)*fes.nDof + fes.ElemsDof;
                Gpert(dof)  = (g_qp{II}.*womegaJ_qp) * psim;
            end
            
        else
            
            Gpert   = G;
            
        end
        
        %Update Gdot:
        Gdot        = Gdot + Udot_scal*(w_deriv(ideriv)*Gpert);
        
    end
    
    %Assemble Gdot + 1/tau*G into F:
    F(model.nDiff*fes.nDof+1:end)   = Gdot + G/tau;
    
    %----------------------------------------------------------------------
    %ASSEMBLY:
    
    if ComputeJ
        
        %Note: there are three blocks: differential variables, algebraic 
        %variables, other variables (such as mesh velocity). Block U
        %contains diff variables, V block contains the rest.
        
        %Jacobian of differential equations:
        F_UV_iv         = cat(1, iv(:), iv_ee(:), iv_ew(:), iv_we(:), iv_ww(:), iv_11(:), iv_22(:));
        F_UV_jv         = cat(1, jv(:), jv_ee(:), jv_ew(:), jv_we(:), jv_ww(:), jv_11(:), jv_22(:));
        F_UV_sv         = cat(1, sv(:), sv_ee(:), sv_ew(:), sv_we(:), sv_ww(:), sv_11(:), sv_22(:));
        F_UV            = sparse(F_UV_iv, F_UV_jv, F_UV_sv, ...
                            model.nDiff*fes.nDof, model.nVars*fes.nDof);
                        
        %Jacobian of restriction:
        G_U                         = sparse(iv_g(:), jv_g(:), sv_g(:), ...
                                        model.nAlg*fes.nDof, model.nDiff*fes.nDof);
        [G_U_iv, G_U_jv, G_U_sv]    = find(G_U);
        
        %H(U,V) := tau*dG/dU(U)*M^{-1}*F(U,V) + G(U)
        %The derivatives of dG/dU are neglected.
        
        %Contribution of dG/dU*M^{-1}*dF/d(U,V):
        M_diag                      = diag(Mm);
        Minv_diag                   = 1.0./M_diag;
        Minv_approx                 = spdiags( repmat(Minv_diag, model.nDiff, 1), 0, ...
                                        model.nDiff*fes.nDof, model.nDiff*fes.nDof );
        Minv_approx                 = MassMatrixExpand(inv(Mm), model.nDiff, 0);
        H_UV                        = G_U * Minv_approx * F_UV;
        [H_UV_iv, H_UV_jv, H_UV_sv] = find(H_UV);
        
        %We need to add the contribution of dG/dU.
        
        %(i,j,s) values of full Jacobian. Here we abuse the notation. F
        %means [F,G], U means [U,V,W]:
        dF_dU.iv    = cat(1, F_UV_iv, ...
                            H_UV_iv+model.nDiff*fes.nDof, ...
                            G_U_iv+model.nDiff*fes.nDof        );
        dF_dU.jv    = cat(1, F_UV_jv, H_UV_jv, G_U_jv );
        dF_dU.sv    = cat(1, F_UV_sv, H_UV_sv, G_U_sv/tau);
        
    else
        dF_dU.iv    = zeros(0,0);
        dF_dU.jv    = zeros(0,0);
        dF_dU.sv    = zeros(0,0);
    end
    
    %----------------------------------------------------------------------
    %Deltat for CFL=1:
    
    %DEBUG:
%     lambdav{1}      = 0.0*lambdav{1};   %Disable mesh distortion
%     lambdav{2}      = 0.0*lambdav{2};   %Disable convection
%     lambdav{3}      = 0.0*lambdav{3};   %Disable diffusion
    
    Deltat_CFL      = Inf;
    for jj=1:length(lambdav)
        %Take Linf norm of lambda at each element:
        lambdav_elems   = max(lambdav{jj}, [], 2);
        Deltat_CFL      = min(Deltat_CFL, min(hp_elems(:).^(jj-1)./lambdav_elems));
    end
    
end

function C=UpsilonMat(A, B)

    m   = size(A,2);
    n   = size(B,2);
    N   = size(A,1);
    C   = zeros(N, m*n);
    for jj=1:n
        for ii=1:m
            C(:, (jj-1)*m+ii)   = A(:,ii).*B(:,jj);
        end
    end
    
end

function [imat, jmat] = J_ElemsDof(ElemsDof1, ElemsDof2)

    n1      = size(ElemsDof1,2);
    n2      = size(ElemsDof2,2);
    nElems      = size(ElemsDof1,1);
    
    imat    = zeros(nElems, n1*n2);
    jmat    = zeros(nElems, n1*n2);
    for ii=1:n1
        for jj=1:n2
            imat(:,ii+(jj-1)*n1)    = ElemsDof1(:,ii);
            jmat(:,ii+(jj-1)*n1)    = ElemsDof2(:,jj);
        end
    end
    
end

    