function uLeg = P1ToQX(unodes, fes)

    uLeg        = zeros(size(fes.ElemsDof));
    u1          = unodes(1:end-1);
    u2          = unodes(2:end);
    uLeg(:,1)   = 0.5*(u1+u2);
    uLeg(:,2)   = 0.5*(u2-u1);
    uLeg        = reshape(uLeg.',fes.mesh.nElems*(fes.p+1),1);
    
end
 