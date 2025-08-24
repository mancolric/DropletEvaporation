%Prolongation matrix from continuous space of degree p to discontinuous
%space of degree P.

function Pm = PXToQX_Matrix(nElems, p)
    
    %Matrix with Legendre polynomials at Chebyshev nodes:
    xiv     = xi_Chebyshev(p);
    phim    = PolyLeg(xiv, p);
    %If p(x) = L_k(x) u_k = phi_j(x) a_j, then
    %   p(x) = phi_j(x_k) L_k(x) a_j 
    %and therefore
    %   u_k = phi_j(x_k) a_j = phim_kj a_j
    %Similarly,
    %   a_i = psi_ij u_j
    %with psi^{-1} = phi:
    smat            = inv(phim);
    [imat,jmat]     = ndgrid(1:p+1, 1:p+1);
    
    %Loop elements and assemble psim:
    Pm_iv           = zeros(p+1, p+1, nElems);
    Pm_jv           = zeros(size(Pm_iv));
    Pm_sv           = zeros(size(Pm_iv));
    for iElem=1:nElems
        %Accumulated dof for elem iElem in QSpace is (iElem-1)*(p+1):
        Pm_iv(:,:,iElem)    = imat + (iElem-1)*(p+1);
        %Accumulated dof for elem iElem in PSpace is (iElem-1)*p:
        Pm_jv(:,:,iElem)    = jmat + (iElem-1)*p;
        %Nnz:
        Pm_sv(:,:,iElem)    = smat;
    end
    
    Pm              = sparse(Pm_iv(:), Pm_jv(:), Pm_sv(:), nElems*(p+1), nElems*p+1);
    
end

    