%Return derivatives (w.r.t. xi) of Legrende polynomials of degree p at points xi.
%The output is a matrix of size (length(xi), p+1):
function dL = dPolyLegendre(xi, p)
    
    %Allocate output:
    xi      = xi(:);
    dL      = zeros(length(xi), p+1);
    L       = PolyLegendre(xi, p);
    
    %Apply recurrence:
    %   dT_(n) = (n)*T_(n-1) + x*dT_(n-1)
    dL(:,1)     = 0.0;
    for nn=1:p
        dL(:,nn+1)     = nn*L(:,nn) + xi.*dL(:,nn);
    end

end
