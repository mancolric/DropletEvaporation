function u1 = ForceP1(u0, fes)

    u1  = P1ToQX(QXToP1(u0,fes), fes);
    
end