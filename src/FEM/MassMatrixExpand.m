function Mm = MassMatrixExpand(M_II, nDiff, nAlg)

    [iv,jv,sv]  = find(M_II);
    [M,N]       = size(M_II);
    
    iv2         = zeros(length(iv),nDiff);
    jv2         = zeros(length(iv),nDiff);
    sv2         = zeros(length(iv),nDiff);
    
    for II=1:nDiff
        iv2(:,II)   = iv + (II-1)*M;
        jv2(:,II)   = jv + (II-1)*N;
        sv2(:,II)   = sv;
    end
    
    Mm          = sparse(iv2(:), jv2(:), sv2(:), (nDiff+nAlg)*M, (nDiff+nAlg)*N);
    
end
