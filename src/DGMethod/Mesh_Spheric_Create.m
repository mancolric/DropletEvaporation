function mesh = Mesh_Spheric_Create(x_faces)

    %Geometry:
    mesh.x_faces    = x_faces(:);       %Always column vector
    mesh.A_faces    = (mesh.x_faces.^2)*(4*pi);
    mesh.V_cells    = ((mesh.x_faces(2:end).^3)-(mesh.x_faces(1:end-1).^3))*(4*pi/3);
    mesh.nElems     = length(mesh.x_faces)-1;
    
    %Jacobian of transformation between physical element and reference
    %element = length of physical element / length of reference element:
    mesh.J          = (mesh.x_faces(2:end)-mesh.x_faces(1:end-1))/2.0;
    
    %Function such that dOmega = omega*dx:
    mesh.omega      = @(r) (4*pi)*(r.^2);
    
end
