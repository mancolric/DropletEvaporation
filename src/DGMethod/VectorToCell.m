function uc=VectorToCell(uv, nVars)

    if ~(isvector(uv) && isnumeric(uv))
        error('First argument must be a vector')
    end
    
    uc      = cell(nVars,1);
    nDof    = length(uv)/nVars;
    for II=1:nVars
        dof     = (II-1)*nDof+1:II*nDof;
        uc{II}  = uv(dof);
    end
    
end
