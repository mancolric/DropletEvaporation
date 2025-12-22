function fes = FES_Create(mesh, p)

    warning('Deprecated. Use FES_QX_Create')
    fes             = FES_QX_Create(mesh, p);
    
end