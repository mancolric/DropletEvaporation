function LUFact = LUFactorization(A)
    
    %Sparse LU factorization:
    %   P*(R\A)*Q = L*U 
    [L,U,P,Q,R]     = lu(A);
    LUFact.L        = L;
    LUFact.U        = U;
    LUFact.P        = P;
    LUFact.Q        = Q;
    LUFact.R        = R;
    
end

    