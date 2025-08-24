function mesh = Mesh_Cartesian_Create(x_faces)

    %Geometry:
    mesh.x_faces    = x_faces(:);       %Always column vector
    mesh.A_faces    = ones(length(x_faces));
    mesh.V_cells    = (mesh.x_faces(2:end)-mesh.x_faces(1:end-1))*mesh.A_faces(1);
    mesh.nElems     = length(mesh.x_faces)-1;
    
    %Jacobian of transformation between physical element and reference
    %element = length of physical element / length of reference element:
    mesh.J          = (mesh.x_faces(2:end)-mesh.x_faces(1:end-1))/2.0;
    
    %Function such that dOmega = omega*dx:
    mesh.omega      = @(x) 0.0*x + mesh.A_faces(1);
    
end
