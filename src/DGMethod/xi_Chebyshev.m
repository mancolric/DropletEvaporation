%Return Chebyshev nodes of degree p at the interval [-1,1]:

function xiv = xi_Chebyshev(p)

    if p==0
        xiv     = [ 0.0 ];
        return
    end
    
    iv      = 0:p;
    xiv     = -cos(pi*iv/p);
    xiv     = xiv(:);
    
end
