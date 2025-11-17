%Prolongation matrix from continuous space of degree p to discontinuous
%space of degree P.

function Pm = PXToQY_Matrix(fesP, fesQ)
    
    %Check number of elements:
    if fesP.mesh.nElems~=fesQ.mesh.nElems
        error('Number of elements associated with both spaces must be the same')
    end
    
    %Load quadrature rules:
    QR      = fesQ.QuadRule;
    
    %u(x) = L_j(x) u_j = phi_j(x) a_j. If we project onto phi_i:
    %   (phi_i, L_j) u_j = (phi_i, phi_j) a_j
    %So 
    %   a = ( phiphi \ phiL ) * u
    phim        = fesQ.NCompute(QR.xi);
    Lm          = fesP.NCompute(QR.xi);
    wm          = diag(QR.w);
    phiL        = phim.'*wm*Lm;
    phiphi      = phim.'*wm*phim;
    prolm_elem  = phiphi\phiL;
    
    %Loop elements and assemble psim:
    nElems          = fesP.mesh.nElems;
    Pm_iv           = zeros(nElems, fesQ.p+1, fesP.p+1);
    Pm_jv           = zeros(size(Pm_iv));
    Pm_sv           = zeros(size(Pm_iv));
    aux_elems       = 1:nElems;
    for iDof=1:fesQ.p+1
        for jDof=1:fesP.p+1
            %Accumulated dof for elem iElem in QSpace is (iElem-1)*(p+1):
            Pm_iv(:,iDof,jDof)  = (aux_elems-1)*(fesQ.p+1) + iDof;
            %Accumulated dof for elem iElem in PSpace is (iElem-1)*p:
            Pm_jv(:,iDof,jDof)  = (aux_elems-1)*fesP.p + jDof;
            %Nnz:
            Pm_sv(:,iDof,jDof)  = prolm_elem(iDof,jDof);
        end
    end
    
    Pm              = sparse(Pm_iv(:), Pm_jv(:), Pm_sv(:), nElems*(fesQ.p+1), nElems*fesP.p+1);
    
end

    