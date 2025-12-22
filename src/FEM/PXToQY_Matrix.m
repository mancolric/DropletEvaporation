%Prolongation matrix from continuous space of degree p to discontinuous
%space of degree P.

function Pm = PXToQY_Matrix(fesP, fesQ)
    
    %Check number of elements:
    if fesP.mesh.nElems~=fesQ.mesh.nElems
        error('Number of elements associated with both spaces must be the same')
    end
    
    %Check spaces orders:
    if fesP.p>fesQ.p
        error('Cannot prolongate onto a space of smaller dimension')
    end
    
    %u(x) = L_j(x) u_j = phi_j(x) a_j. If we evaluate this expresion at the
    %interpolation nodes x_k:
    %   delta_jk u_j = u_k = phi_j(x_k) a_j
    %So 
    %   a = phim(x_nodes) \ u
    prolm_elem                  = zeros(fesQ.p+1, fesP.p+1);
    prolm_elem(1:fesP.p+1,:)    = fesP.Lag_Leg(1:fesP.p+1,:);
    
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

    