%From three pair of data (x_i, y_i), i=1,...,3, obtain parameters 
%   A, C, n from curve y_i = A + C x_i^n.
function [A, C, n] = ExpRegression(xv, yv)

    xv              = xv(:);
    yv              = yv(:);
    [xsorted, I]    = sort(xv);
    ysorted         = yv(I);
    
    %We have
    %   y2-y1                   = C (x2^n-x1^n) = C x1^n ((x2/x1)^n-1)
    %   y3-y2                   = C (x3^n-x2^n) = C x2^n ((x3/x2)^n-1)
    %   ln( (y3-y2)/(y2-y1) )   = n ln(x2/x1) + ln((x3/x2)^n-1) - ln((x2/x1)^n-1) 
    %Initial condition with Richarson extrapolation, assuming x3/x2 =
    %x2/x1:
    n0              = log((ysorted(3)-ysorted(2))/(ysorted(2)-ysorted(1))) / ...
                        log(xsorted(2)/xsorted(1));
    
    function [f,J]  = ResidualFun(n, ComputeJ)
        
        %Residual vector:
        f       = n * log(xsorted(2)/xsorted(1)) + ...
                        log((xsorted(3)/xsorted(2))^n-1) - log((xsorted(2)/xsorted(1))^n-1) - ...
                        log((ysorted(3)-ysorted(2))/(ysorted(2)-ysorted(1)));
        
        %Approximate Jacobian:
        J       = log(xsorted(2)/xsorted(1));
    
    end

    %DEBUG: Plot f(n):
%     nv              = linspace(1.0, 4.0, 100);
%     fv              = zeros(length(nv),1);
%     for ii=1:length(nv)
%         [fv(ii),~]  = ResidualFun(nv(ii), false);
%     end
%     figure()
%     plot(nv, fv)
%     error("aa")

    %Find n with NR method:
    [n, ~, ~]       = NewtonRaphson(@ResidualFun, [n0], 0.0, 1e-3, 100);
    
    %Find other coefficients:
    C               = (ysorted(2)-ysorted(1))/(xsorted(2)^n-xsorted(1)^n);
    A               = ysorted(1) - C*xsorted(1)^n;
    
end

function [A, C, n] = ExpRegression_old2(xv, yv)

    xv  = xv(:);
    yv  = yv(:);
    
    function [f,J]  = ResidualFun(u, ComputeJ)
        
        A   = u(1);
        C   = u(2);
        n   = u(3);
        
        %Residual vector:
        f       = A + C*xv.^n - yv;
        
        %Jacobian:
        J       = zeros(3,3);
        J(:,1)  = 1.0;
        J(:,2)  = xv.^n;
        J(:,3)  = C*n*xv.^(n-1);
        
    end

    %Initial condition: assume yv(imin) is exact:
    [xsorted, I]    = sort(xv);
    ysorted         = yv(I);
    A               = ysorted(1);
    n               = log((ysorted(3)-A)/(ysorted(2)-A)) / ...
                        log(xsorted(3)/xsorted(2));
    C               = (ysorted(3)-A)/xsorted(3)^n;
    
    %Call NR method:
    [usol, ~, ~]    = NewtonRaphson(@ResidualFun, [A; C; n], 0.0, 1e-3, 100);
    A               = usol(1);
    C               = usol(2);
    n               = usol(3);
    
end
   
function [A, C, n] = ExpRegression_old(xv, yv)

    [xsorted, I]    = sort(xv);
    ysorted         = yv(I);
    
    %If y = A + C*x^n, and x is small, between two experiments we will have
    %   Delta y ~= C n (Delta x)^(n-1) + O((Delta x)^n)
    %Therefore:
    n               = 1 + log((ysorted(3)-ysorted(2))/(ysorted(2)-ysorted(1))) / ...
                            log((xsorted(3)-xsorted(2))/(xsorted(2)-xsorted(1)));
    C               = (ysorted(2)-ysorted(1)) / (  n*(xsorted(2)-xsorted(1))^(n-1)   );
    A               = ysorted(1) - C*xsorted(1)^n;

end