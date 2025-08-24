function P = Polyfit(x,y,n)

    [p, ~, mu]  = polyfit(x,y,n);
    
    P.p         = p;
    P.mu        = mu;
    
end