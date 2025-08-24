%   Sm = ScalMatrix(sv, N)
%Create a block-diagonal scaling matrix. The diagonal factors for each
%block are stored in sv and the size of each block is N:
function Sm = ScalMatrix(sv, N)

    nBlocks     = length(sv);
    scalv       = zeros(nBlocks*N,1);
    for II=1:nBlocks
        scalv((II-1)*N+1:II*N)  = sv(II);
    end
    Sm          = spdiags(scalv, [0], nBlocks*N, nBlocks*N);
    
end    