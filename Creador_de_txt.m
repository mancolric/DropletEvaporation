% Creador_de_txt.m
% Script para generar el archivo Parametros.txt
% Define aquí las variables, y el script generará el txt automáticamente.

clear; clc;

% 1. Definición de Parámetros

% Identificador de Simulación
Id_sim            = 1009;

% Discretización
nElems_l          = 12;          % Number of elements in the liquid phase
nElems_g          = 500;         % Number of elements in the gas phase
hmin_l            = 1.00E-06;
hmin_g            = 1.00E-09;
p                 = 5;           % Degree of the polynomial basis functions

Deltat0           = 1.00E-09;    % Initial time step
t_final           = 1.00E+02;    % Final simulation time
TimeAdapt         = 1;           % Enable adaptive time-stepping (1 = true, 0 = false)
TolT              = 1.00E-04;    % Tolerance for temporal error

% Guardados
n_saved_solutions = 2;           % Number of solution snapshots to save throughout simulation
n_saves           = 2;           % Number of intermediate "safety saves" during simulation

% Condiciones Iniciales y CC
R_0               = 5.00E-04;    % Initial droplet radius [m]
XRad              = 150;         % Domain size (in gas phase) as a multiple of the initial radius
R_end_percent     = 0.2;         % Final droplet size as a fraction of initial radius

T_0               = 300;         % Initial temperature at the center of the droplet [K]
T_inf             = 1600;        % Ambient (far-field) temperature [K]

% Especies
fuel_names        = 'Heptano, Etanol';
mass_fracL        = '0.5, 0.5';
inert_comps       = 'N2,O2,CO2,H2O';
mass_fracG        = '0.79,0.21,0.0,0.0';


% COMBUSTIÓN
PreExp            = '6.3e11, 3.981e14, 5e8';
ActEnergy         = '125520, 167360, 167360';
T_exp             = '0, 0, 0';

DiffStoi          = [0, -11, 7, 8, -1, 0; % N, O, CO2, H2O, Fuels
    0, -3, 2, 3, 0, -1];


ROrder            = [0, 1, 0, 0, 1, 0; % N, O, CO2, H2O, Fuels
    0, 1, 0, 0, 0, 1];


%2. Generación del Archivo Parametros.txt

nombre_archivo = 'Parametros.txt';
fileID = fopen(nombre_archivo, 'w');

if fileID == -1
    error('No se pudo crear el archivo %s.', nombre_archivo);
end

% Agrupamos todos los datos en el orden requerido por Run_Cluster.m
% mat2str convierte las matrices a texto para que se guarden correctamente en el txt
datos = {Id_sim, nElems_l, nElems_g, hmin_l, hmin_g, p, Deltat0, t_final, TimeAdapt, TolT, ...
    n_saved_solutions, n_saves, R_0, XRad, R_end_percent, T_0, T_inf, fuel_names, ...
    mass_fracL, inert_comps, mass_fracG, PreExp, ActEnergy, T_exp, mat2str(DiffStoi), mat2str(ROrder)};

num_columnas = length(datos);

for j = 1:num_columnas
    valor = datos{j};
    
    if ischar(valor) || isstring(valor)
        fprintf(fileID, '%s', valor);
    elseif isnumeric(valor)
        % Si es entero, se imprime como entero
        if mod(valor, 1) == 0 && abs(valor) < 1e6
            fprintf(fileID, '%d', valor);
        else
            % Columnas que normalmente llevan formato científico:
            % 4:hmin_l, 5:hmin_g, 7:Deltat0, 8:t_final, 10:TolT, 13:R_0
            if ismember(j, [4, 5, 7, 8, 10, 13])
                fprintf(fileID, '%.2E', valor);
            else
                fprintf(fileID, '%g', valor);
            end
        end
    end
    
    % Añadir un tabulador entre cada columna, excepto en la última
    if j < num_columnas
        fprintf(fileID, '\t');
    end
end

% Finalizar con un salto de línea
fprintf(fileID, '\n');
fclose(fileID);

fprintf('¡Éxito! El archivo "%s" ha sido generado.\n', nombre_archivo);
