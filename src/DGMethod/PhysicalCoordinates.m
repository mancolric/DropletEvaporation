%Given some coordinates xi in the reference element, return a matrix 
%(mesh.Nx, length(xi)) with the corresponding coordinates at the physical elements.
function x = PhysicalCoordinates(mesh, xi)

    %For a given element "i", the quadrature nodes are placed at
    % x     = x_faces(i)*0.5*(1-xi) + x_faces(i+1)*0.5*(1+xi). 
    xi      = xi(:);    %Always column
    x       = mesh.x_faces(1:end-1)*0.5*(1.0-xi') + ...
                mesh.x_faces(2:end)*0.5*(1.0+xi');
            
end