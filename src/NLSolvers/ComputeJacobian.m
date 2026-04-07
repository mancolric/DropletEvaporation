%Compute Jacobian of function f using centered finite differences. Use only
%with small matrices for checking purposes.
function J=ComputeJacobian(fun, u0, delta)

Delta_vector = u0*delta + delta/2;
    %Get size of f:
    f0          = fun(u0);
    Nf          = length(f0);
    
    %Compute Jacobian with finite differences:
    J           = zeros(Nf);
    for jj=1:Nf
        %Compute f(u-delta*e_j)
        u1      = u0;
        u1(jj)  = u1(jj)-Delta_vector(jj);
        f1      = fun(u1);
        %Compute f(u+delta*e_j):
        u2      = u0;
        u2(jj)  = u2(jj)+Delta_vector(jj);
        f2      = fun(u2);
        %Apply centered finite differences:
        J(:,jj) = (f2-f1)/(2*Delta_vector(jj));
    end
    
end