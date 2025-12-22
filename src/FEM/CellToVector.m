function uv = CellToVector(uc)

    if ~iscell(uc)
        error('Argument must be a cell')
    end
    
    nVars   = length(uc);
    nDof    = length(uc{1});
    uv      = zeros(nVars*nDof,1);
    for II=1:nVars
        dof         = (II-1)*nDof+1:II*nDof;
        uv(dof)     = uc{II};
    end
    
end
