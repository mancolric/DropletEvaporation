%[c1, c2, ...] = Cells_Allocate(m, n, ComputeJ, u)
%
%Returns ci=cell(m,n). If ComputeJ==true, ci{ii,jj}=zeros(size(u));
%otherwise, the cell is empty.

function varargout = Cells_Allocate(m, n, ComputeJ, u)

    for II=1:nargout
        var_cell        = cell(m,n);
        if ComputeJ
            for ii=1:m
                for jj=1:n
                    var_cell{ii,jj} = zeros(size(u));
                end
            end
        end
        varargout{II}   = var_cell;
    end
    
end
