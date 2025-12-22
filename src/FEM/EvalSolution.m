%Given some coordinates xi in the reference element, return a matrix 
%(mesh.Nx, length(xi)) with the corresponding values of the solution:
function u = EvalSolution(usol, fes, elems, xi)

    nVars   = length(usol);
    
    %Matrix (N, p+1) with the values of the shape functions:
    L       = fes.NCompute(xi);
    
    %Evaluate solution at natural coordinates xi:
    u           = cell(nVars,1);
    ElemsDof    = fes.ElemsDof(elems,:);
    for II=1:nVars
        %Matrix (p+1, Nx) with the solution at each element:
        u_n_mat = reshape( usol{II}(ElemsDof), size(ElemsDof) ); %reshape necessary in case nElems=1
        %Values of the solution at each element and each value of xi:
        u{II}   = u_n_mat*L.';
    end
    
end