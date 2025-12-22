%Return Legrende polynomials of degree p at points xi. The output is a
%matrix of size (length(xi), p+1):
function L = PolyLegendre(xi, p)
    
    %Allocate output:
    xi      = xi(:);
    L       = zeros(length(xi), p+1);
    
    %Apply recurrence:
    %   T0      = 1
    %   T1      = x
    %   T_(n+1) = ((2n+1)x*T_n - n*T_(n-1))/(n+1)
    L(:,1)      = 1.0;
    if p>0
        L(:,2)  = xi;
        for nn=1:p-1
            L(:,nn+2)   = ((2.0*nn+1).*xi.*L(:,nn+1) - nn*L(:,nn))/(nn+1);
        end
    end
    
end

