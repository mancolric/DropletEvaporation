function Run_Cluster(Id_sim)
format long            % Display numerical results with high precision
close all              % Close all Plots
clc                    % Clear Command Window

ruta_actual       = fileparts(mfilename('fullpath'));
data = readcell(fullfile(ruta_actual,'Parametros.txt'), 'Delimiter', '\t');% MAIN SCRIPT TO LAUNCH THE DROPLET EVAPORATION SIMULATION
% data = readmatrix(fullfile('/home/tfd/dbetran/Simulaciones_14000','Parametros.txt'));
simID             = [data{Id_sim, 1}];


addpath(genpath(pwd)); % Add all subfolders of the current directory to the MATLAB path

% ---------------------- MESH AND TIME CONFIGURATION ----------------------

nElems_l          = [data{Id_sim, 2}];    % Number of elements in the liquid phase
nElems_g          = [data{Id_sim, 3}];    % Number of elements in the gas phase
hmin_l            = [data{Id_sim, 4}];  
hmin_g            = [data{Id_sim, 5}]; 
p                 = [data{Id_sim, 6}];     % Degree of the polynomial basis functions

Deltat0           = [data{Id_sim, 7}];  % Initial time step
                           % If you get warnings like:
                           %   “Nonlinear solver did not converge at stage 2. Reducing time step”
                           %   or
                           %   "Time error (0.xxxxx) too large. Multiplying time step by a factor 0.xxxxx",
                           % you may reduce Deltat0 to avoid them
                           % However, this is optional — time step control will automatically handle it

t_final           = [data{Id_sim, 8}];   % Final simulation time
                           % Use the actual time you want the simulation to stop at
                           % If you want the simulation to continue until the droplet disappears, 
                           % set a value larger than the expected final time

TimeAdapt         = [data{Id_sim, 9}];  % Enable adaptive time-stepping (true/false)
TolT              = [data{Id_sim, 10}];  % Tolerance for temporal error (<1e-3) (used if TimeAdapt = true)

% ---------------------------- OUTPUT OPTIONS -----------------------------

PlotRes           = true;  % Plot intermediate results (true/false)
                           % NOTE: if true, Save must be false

Save              = false; % Save results to file (true/false)
                           % NOTE: if true, PlotRes must be false
                           
n_saved_solutions = [data{Id_sim, 11}];   % Number of solution snapshots to save throughout simulation
n_saves           = [data{Id_sim, 12}];    % Number of intermediate "safety saves" during simulation

% ---------------------- INITIAL AND BOUNDARY CONDITIONS ------------------

R_0               = [data{Id_sim, 13}];    % Initial droplet radius [m]
XRad              = [data{Id_sim, 14}];        % Domain size (in gas phase) as a multiple of the initial radius
R_end_percent     = [data{Id_sim, 15}];       % Final droplet size as a fraction of initial radius (when Save = true)

T_0               = [data{Id_sim, 16}];        % Initial temperature at the center of the droplet [K]
T_inf             = [data{Id_sim, 17}];       % Ambient (far-field) temperature [K]

% --------------------------- SPECIES DEFINITION --------------------------

% fuel_names        = %{'Heptano'};              % Fuel components
fuel_names        = strsplit([data{Id_sim, 18}], ',');
mass_fracL        = num2cell(str2double(strsplit(string(data{Id_sim, 19}), ',')));
                                                            % Mass fraction of each fuel component (same order as above)
inert_comps       = {'N2', 'O2', 'CO2', 'H2O'};             % Inert gas species (DO NOT MODIFY)

% DE AQUÍ SOLO SE USA EL % DE N2 EN COMBUSTIÓN
mass_fracG        = num2cell(str2double(strsplit([data{Id_sim, 20}], ',')));
                                                            % Mass fraction of each inert species in the gas phase (same order as above)

% --------------------------- COMBUSTION PARAMET --------------------------

PreExp            = [data{Id_sim, 21}];
ActEnergy         = [data{Id_sim, 22}];
T_exp             = [data{Id_sim, 23}];
Fuel_exp          = [data{Id_sim, 24}];
O2_exp            = [data{Id_sim, 25}];

% ------------------------ INITIAL SIMULATION CALL ------------------------

% Run the full droplet evaporation simulation from t = 0
folderName  = sprintf('Sim_%d', simID);
mkdir(ruta_actual, folderName);
[~] = droplet_test_Cluster(nElems_l, nElems_g, hmin_l, hmin_g, p, Deltat0, t_final, ...
    TimeAdapt, TolT, PlotRes, Save, fuel_names, mass_fracL, inert_comps, ...
    mass_fracG, n_saved_solutions, T_0, T_inf, R_0, XRad, R_end_percent, ...
    n_saves, folderName, PreExp, ActEnergy, T_exp, Fuel_exp, O2_exp);
return

end