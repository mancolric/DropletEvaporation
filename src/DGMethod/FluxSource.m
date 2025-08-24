function [b, J] = FluxSource(model, t, usol, fes, ComputeJ)

    mesh        = fes.mesh;
    p           = fes.p;
    QuadRule    = fes.QuadRule;
    nElems      = mesh.nElems;
    
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
    nVars       = length(u_qp);
    
    %Evaluate flux and source terms at quadrature points:
    [f_qp, df_du_qp, df_dgradu_qp, ...
        Q_qp, dQ_du_qp, dQ_dgradu_qp]   = model.fQ(model, t, x_qp, u_qp, du_dx_qp, ComputeJ);
    
    %Load shape functions in reference element at quadrature nodes:
    %Matrix is of size (#quadrature nodes, p+1):
    psim        = PolyLeg(QuadRule.xi, p);
    dpsim       = dPolyLeg(QuadRule.xi, p);
    
    %Allocate rhs vector:
    % b_(I,alpha)   = int(Q_I*psi_alpha) + int(f_I*dpsi_alpha/dx) - bint(fhat_I*psi_alpha)
    b           = zeros(nVars*fes.nDof, 1);
    
    %Contribution of flux and source terms in the domain:
    for II=1:nVars
        b((II-1)*fes.nDof+fes.ElemsDof)     = (f_qp{II}.*Jinv_qp.*womegaJ_qp) * dpsim + ...
                                                (Q_qp{II}.*womegaJ_qp) * psim;
    end
    
    %Compute db_(I,alpha)/du_(J,beta). 
    if ComputeJ
        
        %Allocate memory for jacobian elements. For each (II,JJ) block, we have
        %nElems*(p+1)^2 nonzeros:
        iv      = zeros(mesh.nElems, (p+1)*(p+1), nVars, nVars);
        jv      = zeros(size(iv));
        sv      = zeros(size(iv));
        
        %Products of shape functions:
        UpsilonNN       = UpsilonMat(psim, psim);
        UpsilonGN       = UpsilonMat(dpsim, psim);
        UpsilonNG       = UpsilonMat(psim, dpsim);
        UpsilonGG       = UpsilonMat(dpsim, dpsim);
       
        %Loop blocks:
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof, fes.ElemsDof);
        for II=1:nVars
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
    %FLUXES AT INTERNAL FACES:
    
    %Solution at west and east faces:
    u_w         = EvalSolution(usol, fes, 2:nElems, [-1.0]); %cell(nVars,1)
    du_dx_w     = EvalSolution_dx(usol, fes, 2:nElems, [-1.0]); %cell(nVars,1)
    u_e         = EvalSolution(usol, fes, 1:nElems-1, [1.0]); %cell(nVars,1)
    du_dx_e     = EvalSolution_dx(usol, fes, 1:nElems-1, [1.0]); %cell(nVars,1)
    psim_w      = PolyLeg([-1.0], p);
    psim_e      = PolyLeg([1.0], p);
    dpsim_w     = dPolyLeg([-1.0], p);
    dpsim_e     = dPolyLeg([1.0], p);
    Jinv_elems  = 1.0./mesh.J;
    omega_faces = mesh.omega(mesh.x_faces(2:end-1));
    h_elems     = diff(mesh.x_faces);
    h_faces     = 2.0./(1.0./h_elems(1:nElems-1)+1.0./h_elems(2:nElems));
    
    %Fluxes at the nElems-1 internal faces:
    [fhat, dfhat_du_e, dfhat_dgradu_e, ...
        dfhat_du_w, dfhat_dgradu_w]     = model.ftilde(model, t, mesh.x_faces(2:end-1), ...
                                            u_e, du_dx_e, u_w, du_dx_w, h_faces, ...
                                            ComputeJ);
                                        
    %Fluxes are leaving the elements 1:nElems-1 at the east faces:
    for II=1:nVars
        dof         = (II-1)*fes.nDof+fes.ElemsDof(1:nElems-1,:);
        b(dof)      = b(dof) - (fhat{II}.*omega_faces) * psim_e;
    end
    if ComputeJ
        
        %Variations of leaving fluxes w.r.t. variations at left:
        iv_ee           = zeros(nElems-1, (p+1)*(p+1), nVars, nVars);
        jv_ee           = zeros(size(iv_ee));
        sv_ee           = zeros(size(iv_ee));
        UpsilonNN       = UpsilonMat(psim_e, psim_e);
        UpsilonNG       = UpsilonMat(psim_e, dpsim_e);
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof(1:nElems-1,:), fes.ElemsDof(1:nElems-1,:));
        for II=1:nVars
            for JJ=1:nVars
                iv_ee(:,:,II,JJ)    = imat + (II-1)*fes.nDof;
                jv_ee(:,:,II,JJ)    = jmat + (JJ-1)*fes.nDof;
                sv_ee(:,:,II,JJ)    = - (dfhat_du_e{II,JJ}.*omega_faces) * UpsilonNN - ...
                                        (dfhat_dgradu_e{II,JJ}.*Jinv_elems(1:nElems-1).*omega_faces) * UpsilonNG;
            end
        end
       
        %Variations of leaving fluxes w.r.t. variations at the right:
        iv_ew           = zeros(nElems-1, (p+1)*(p+1), nVars, nVars);
        jv_ew           = zeros(size(iv_ew));
        sv_ew           = zeros(size(iv_ew));
        UpsilonNN       = UpsilonMat(psim_e, psim_w);
        UpsilonNG       = UpsilonMat(psim_e, dpsim_w);
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof(1:nElems-1,:), fes.ElemsDof(2:nElems,:));
        for II=1:nVars
            for JJ=1:nVars
                iv_ew(:,:,II,JJ)    = imat + (II-1)*fes.nDof;
                jv_ew(:,:,II,JJ)    = jmat + (JJ-1)*fes.nDof;
                sv_ew(:,:,II,JJ)    = - (dfhat_du_w{II,JJ}.*omega_faces) * UpsilonNN - ...
                                        (dfhat_dgradu_w{II,JJ}.*Jinv_elems(2:nElems).*omega_faces) * UpsilonNG;
            end
        end
        
    end
    
    %Fluxes are entering the elements 2:nElems at the west faces:
    for II=1:nVars
        dof         = (II-1)*fes.nDof+fes.ElemsDof(2:nElems,:);
        b(dof)      = b(dof) + (fhat{II}.*omega_faces) * psim_w;
    end
    if ComputeJ
        
        %Variations of entering fluxes w.r.t. variations at left:
        iv_we           = zeros(nElems-1, (p+1)*(p+1), nVars, nVars);
        jv_we           = zeros(size(iv_we));
        sv_we           = zeros(size(iv_we));
        UpsilonNN       = UpsilonMat(psim_w, psim_e);
        UpsilonNG       = UpsilonMat(psim_w, dpsim_e);
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof(2:nElems,:), fes.ElemsDof(1:nElems-1,:));
        for II=1:nVars
            for JJ=1:nVars
                iv_we(:,:,II,JJ)    = imat + (II-1)*fes.nDof;
                jv_we(:,:,II,JJ)    = jmat + (JJ-1)*fes.nDof;
                sv_we(:,:,II,JJ)    = + (dfhat_du_e{II,JJ}.*omega_faces) * UpsilonNN + ...
                                        (dfhat_dgradu_e{II,JJ}.*Jinv_elems(1:nElems-1).*omega_faces) * UpsilonNG;
            end
        end
       
        %Variations of entering fluxes w.r.t. variations at the right:
        iv_ww           = zeros(nElems-1, (p+1)*(p+1), nVars, nVars);
        jv_ww           = zeros(size(iv_ww));
        sv_ww           = zeros(size(iv_ww));
        UpsilonNN       = UpsilonMat(psim_w, psim_w);
        UpsilonNG       = UpsilonMat(psim_w, dpsim_w);
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof(2:nElems,:), fes.ElemsDof(2:nElems,:));
        for II=1:nVars
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
    [f1, df1_du1, df1_du1_dx] = model.ftilde1(model, t, mesh.x_faces(1), ...
                                    u1, du1_dx, h_elems(1), ComputeJ);
    for II=1:nVars
        dof     = fes.ElemsDof(1,:) + (II-1)*fes.nDof;
        b(dof)  = b(dof) + reshape((f1{II}.*omega_face1) * psim_w, [], 1); %reshape is necessary because it is only one element
    end
    if ComputeJ
        %Variations of incoming flux at face 1 w.r.t. variations in
        %solution at face 1:
        iv_11           = zeros(1, (p+1)*(p+1), nVars, nVars);
        jv_11           = zeros(size(iv_11));
        sv_11           = zeros(size(iv_11));
        UpsilonNN       = UpsilonMat(psim_w, psim_w);
        UpsilonNG       = UpsilonMat(psim_w, dpsim_w);
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof(1,:), fes.ElemsDof(1,:));
        for II=1:nVars
            for JJ=1:nVars
                iv_11(:,:,II,JJ)    = imat + (II-1)*fes.nDof;
                jv_11(:,:,II,JJ)    = jmat + (JJ-1)*fes.nDof;
                sv_11(:,:,II,JJ)    = (df1_du1{II,JJ} .* omega_face1) * ...
                                        UpsilonNN + ...
                                        (df1_du1_dx{II,JJ} .* Jinv_elems(1) .* omega_face1) * ...
                                        UpsilonNG;
            end
        end
                
    end
    
    %Right boundary:
    uN              = EvalSolution(usol, fes, [nElems], [1.0]); %cell(nVars,1)
    duN_dx          = EvalSolution_dx(usol, fes, [nElems], [1.0]); %cell(nVars,1)
    omega_faceN     = mesh.omega(mesh.x_faces(nElems+1));
    
    %Flux at right boundary is leaving the last element:
    [fN, dfN_duN, dfN_duN_dx] = model.ftildeN(model, t, mesh.x_faces(end), ...
                                    uN, duN_dx, h_elems(nElems), ComputeJ);
    for II=1:nVars
        dof     = fes.ElemsDof(nElems,:) + (II-1)*fes.nDof;
        b(dof)  = b(dof) - reshape((fN{II}.*omega_faceN) * psim_e, [], 1); %reshape is necessary because it is only one element
    end
    if ComputeJ
        %Variations of incoming flux at face N w.r.t. variations in
        %solution at face N:
        iv_22           = zeros(1, (p+1)*(p+1), nVars, nVars);
        jv_22           = zeros(size(iv_22));
        sv_22           = zeros(size(iv_22));
        UpsilonNN       = UpsilonMat(psim_e, psim_e);
        UpsilonNG       = UpsilonMat(psim_e, dpsim_e);
        [imat, jmat]    = J_ElemsDof(fes.ElemsDof(nElems,:), fes.ElemsDof(nElems,:));
        for II=1:nVars
            for JJ=1:nVars
                iv_22(:,:,II,JJ)    = imat + (II-1)*fes.nDof;
                jv_22(:,:,II,JJ)    = jmat + (JJ-1)*fes.nDof;
                sv_22(:,:,II,JJ)    = - (dfN_duN{II,JJ} .* omega_faceN) * ...
                                        UpsilonNN - ...
                                        (dfN_duN_dx{II,JJ} .* Jinv_elems(nElems) .* omega_faceN) * ...
                                        UpsilonNG;
            end
        end
                
    end
    
    if ComputeJ
        J       = sparse(   cat(1, iv(:), iv_ee(:), iv_ew(:), iv_we(:), iv_ww(:), iv_11(:), iv_22(:)), ...
                            cat(1, jv(:), jv_ee(:), jv_ew(:), jv_we(:), jv_ww(:), jv_11(:), jv_22(:)), ...
                            cat(1, sv(:), sv_ee(:), sv_ew(:), sv_we(:), sv_ww(:), sv_11(:), sv_22(:)) );
    else
        J       = zeros(0,0);
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

    