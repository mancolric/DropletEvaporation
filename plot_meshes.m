load("mesh_lineal.mat")
figure(1)
plot(xmesh_g, zeros(size(xmesh_g)), 'r*', 'MarkerSize', 6, 'LineWidth', 1.5); 
hold on
load("mesh_exponencial_115.mat")
plot(xmesh_g, zeros(size(xmesh_g)), 'b.', 'MarkerSize', 6, 'LineWidth', 1.5);

legend("Mesh Lineal", "Mesh Combustion")
xlim([0.001068-0.01*0.0010682, 0.001068+0.01*0.0010682])