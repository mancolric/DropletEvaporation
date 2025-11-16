function x=LUSolve(F, b)

    %   A x                     = b
    %   P R^{-1} A Q Q^{-1} x   = P R^{-1} b
    %   L U Q^{-1} x            = P R^{-1} b
    %   x                       = Q U^{-1} L^{-1} * (P R^{-1} b)
    x   = F.Q * ( F.U \ (F.L \ ( F.P * (F.R.\b) ) ) );
    
end