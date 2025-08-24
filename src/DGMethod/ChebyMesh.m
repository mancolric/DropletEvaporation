%Create Chebyshev mesh:

function xnodes = ChebyMesh(xverts, p)

    nElems      = length(xverts)-1;
    
    xnodes      = zeros(nElems*p+1,1);
    xiv         = xi_Chebyshev(p);
    xiv_inner   = xiv(2:p);
    for iElem=1:nElems
        nodes           = (iElem-1)*p+2:(iElem-1)*p+p;
        xnodes(nodes)   = 0.5*(xverts(iElem)+xverts(iElem+1)) + ...
                            0.5*(xverts(iElem+1)-xverts(iElem))*xiv_inner;
    end
    xnodes(1:p:nElems*p+1)  = xverts(:);
    
end