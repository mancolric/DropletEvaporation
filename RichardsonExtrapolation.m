function RichardsonExtrapolation(Tv, r)

    alpha   = log((Tv(3)-Tv(2))/(Tv(2)-Tv(1)))/log(r);
    C       = (Tv(2)-Tv(1))/(r^alpha-1);
    T0      = Tv(1)-C;
    
    figure()
    hplot   = linspace(1,r^2,100);
    Tplot   = T0 + C*hplot.^alpha;
    loglog(hplot, Tplot-T0, 'b')
    hold on
    loglog([1,r,r^2], Tv-T0, 'ob')
    display(alpha)
    display(T0)
    
end