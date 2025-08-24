function [F, dF_dU, dF_dq1, dF_dqN] = CG_fgQ(model, t, usol, fes, ComputeJ)

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
        g_qp, dg_du_qp, dg_dgradu_qp]   = model.fQg(model, t, x_qp, u_qp, du_dx_qp, ComputeJ);
    
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
    
%     if ComputeJ
%         J       = sparse(iv(:), jv(:), sv(:));
%     else
%         J       = zeros(0,0);
%     end
%     
%     return
    
    %----------------------------------------------------------------------
    %FLUXES AT BOUNDARY FACES (I.E., IMPOSITION OF B.C.):
    
    %Left boundary:
    u1              = EvalSolution(usol, fes, [1], [-1.0]); %cell(nVars,1)
    du1_dx          = EvalSolution_dx(usol, fes, [1], [-1.0]); %cell(nVars,1)
    omega_face1     = mesh.omega(mesh.x_faces(1));
    psim_w          = fes.NCompute([-1.0]);
    dpsim_w         = fes.dNCompute([-1.0]);
    
    %Flux at left boundary is entering the first element:
    [f1, df1_du1, df1_du1_dx, df1_dq1] = ...
        model.ftilde1(model, t, mesh.x_faces(1), ...
                        u1, du1_dx, (mesh.x_faces(2)-mesh.x_faces(1))/p, ComputeJ);
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
                                        (df1_du1_dx{II,JJ} ./ mesh.J(1) .* omega_face1) * ...
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
    psim_e          = fes.NCompute([1.0]);
    dpsim_e         = fes.dNCompute([1.0]);
    
    %Flux at right boundary is leaving the last element:
    [fN, dfN_duN, dfN_duN_dx, dfN_dqN] = ...
        model.ftildeN(model, t, mesh.x_faces(end), ...
                        uN, duN_dx, (mesh.x_faces(end)-mesh.x_faces(end-1))/p, ComputeJ);
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
                                        (dfN_duN_dx{II,JJ} ./ mesh.J(end) .* omega_faceN) * ...
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
    
    %int( g psi ):
    for II=1:model.nAlg
        dof     = (model.nDiff+II-1)*fes.nDof + fes.ElemsDof;
        F(dof)  = (g_qp{II}.*womegaJ_qp) * psim;
    end
    if ComputeJ
        
        %Allocate memory:
        iv_g    = zeros(mesh.nElems, (p+1)*(p+1), model.nAlg, nVars);
        jv_g    = zeros(size(iv_g));
        sv_g    = zeros(size(iv_g));
        
        %Products of shape functions:
        UpsilonNN       = UpsilonMat(psim, psim);
        UpsilonNG       = UpsilonMat(psim, dpsim);
        
        %int(dg/du psi psi + dg/d(du/dx) dpsi/dx psi):
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof, fes.ElemsDof);
        for II=1:model.nAlg
            for JJ=1:nVars
                iv_g(:,:,II,JJ)     = imat + (model.nDiff+II-1)*fes.nDof;
                jv_g(:,:,II,JJ)     = jmat + (JJ-1)*fes.nDof;
                sv_g(:,:,II,JJ)     = (dg_du_qp{II,JJ}.*womegaJ_qp) * UpsilonNN + ...
                                        (dg_dgradu_qp{II,JJ}.*Jinv_qp.*womegaJ_qp) * UpsilonNG;
            end
        end
        
    end
    
    %----------------------------------------------------------------------
    %ASSEMBLY:
    
    if ComputeJ
        dF_dU.iv    = cat(1, iv(:), iv_11(:), iv_22(:), iv_g(:));
        dF_dU.jv    = cat(1, jv(:), jv_11(:), jv_22(:), jv_g(:));
        dF_dU.sv    = cat(1, sv(:), sv_11(:), sv_22(:), sv_g(:));                    
    else
        dF_dU.iv    = zeros(0,0);
        dF_dU.jv    = zeros(0,0);
        dF_dU.sv    = zeros(0,0);
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

    