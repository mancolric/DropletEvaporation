%Given some coordinates xi in the reference element, return a matrix 
%(mesh.Nx, length(xi)) with the corresponding values of the derivatives of 
%solution w.r.t. the physical variable "x":
function du_dx = EvalSolution_dx(usol, fes, elems, xi)

    p       = fes.p;
    nVars   = length(usol);
    
    %Matrix (N, p+1) with the values of the derivatives of the shape functions
    %w.r.t. xi:
    dL_dxi  = fes.dNCompute(xi);
    
    %Derivatives w.r.t. x:
    du_dx       = cell(nVars,1);
    ElemsDof    = fes.ElemsDof(elems,:);
    for II=1:nVars
        du_dx{II}       = zeros(length(elems), length(xi));
        u_n_mat         = reshape( usol{II}(ElemsDof), size(ElemsDof) ); %reshape necessary in case nElems=1
        du_dxi          = u_n_mat*dL_dxi.';
        for jj=1:length(xi)
            du_dx{II}(:,jj)     = du_dxi(:,jj)./fes.mesh.J(elems);
        end
    end
    
end
