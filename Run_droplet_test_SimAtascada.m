% MAIN SCRIPT TO LAUNCH THE DROPLET EVAPORATION SIMULATION

format long            % Display numerical results with high precision
clear                  % Clear all variables from workspace
% clc                    % Clear command window
% close all              % Close all figure windows

addpath(genpath(pwd)); % Add all subfolders of the current directory to the MATLAB path

% ---------------------- MESH AND TIME CONFIGURATION ----------------------

nElems_l          = 10;    % Number of elements in the liquid phase
nElems_g          = 12;    % Number of elements in the gas phase
hmin_l            = 1e-6;  
hmin_g            = 1e-6; 
p                 = 5;     % Degree of the polynomial basis functions

Deltat0           = 1e-6;  % Initial time step
                           % If you get warnings like:
                           %   “Nonlinear solver did not converge at stage 2. Reducing time step”
                           %   or
                           %   "Time error (0.xxxxx) too large. Multiplying time step by a factor 0.xxxxx",
                           % you may reduce Deltat0 to avoid them
                           % However, this is optional — time step control will automatically handle it

t_final           = 1e2;   % Final simulation time
                           % Use the actual time you want the simulation to stop at
                           % If you want the simulation to continue until the droplet disappears, 
                           % set a value larger than the expected final time

TimeAdapt         = true;  % Enable adaptive time-stepping (true/false)
TolT              = 1e-4;  % Tolerance for temporal error (<1e-3) (used if TimeAdapt = true)

% ---------------------------- OUTPUT OPTIONS -----------------------------

PlotRes           = false;  % Plot intermediate results (true/false)
                           % NOTE: if true, Save must be false

Save              = true; % Save results to file (true/false)
                           % NOTE: if true, PlotRes must be false
                           
n_saved_solutions = 500;   % Number of solution snapshots to save throughout simulation
n_saves           = 2;    % Number of intermediate "safety saves" during simulation

% ---------------------- INITIAL AND BOUNDARY CONDITIONS ------------------

R_0               = 0.00025; % Initial droplet radius [m]
XRad              = 150;        % Domain size (in gas phase) as a multiple of the initial radius
R_end_percent     = 0.2;       % Final droplet size as a fraction of initial radius (when Save = true)

T_0               = 300;        % Initial temperature at the center of the droplet [K]
T_inf             = 1730;       % Ambient (far-field) temperature [K]

% --------------------------- SPECIES DEFINITION --------------------------

fuel_names        = {'Heptano', 'Butanol'};                    % Fuel components
mass_fracL        = {0.6, 0.4};                         % Mass fraction of each fuel component (same order as above)
inert_comps       = {'N2', 'O2', 'CO2', 'H2O'};     % Inert gas species (DO NOT MODIFY)
mass_fracG        = {0.79, 0.0, 0.0, 0.21};  % Mass fraction of each inert species in the gas phase (same order as above)

% ------------------------ INITIAL SIMULATION CALL ------------------------

% Run the full droplet evaporation simulation from t = 0
[save_vars] = droplet_test(nElems_l, nElems_g, hmin_l, hmin_g, p, Deltat0, t_final, ...
    TimeAdapt, TolT, PlotRes, Save, fuel_names, mass_fracL, inert_comps, ...
    mass_fracG, n_saved_solutions, T_0, T_inf, R_0, XRad, R_end_percent, ...
    n_saves);
return

% --------------------------- RESTART SIMULATION --------------------------

% NOTE: Deltat0 may need to be smaller
% If the simulation is performed in parts, this call restarts it from 
% a previously saved state (e.g., after an interruption or for long simulations)
FileName    = 'Saved_solutions4_51.mat'; % File containing saved variables from previous run

[save_vars] = droplet_testRestart(nElems_l, nElems_g, p, Deltat0, ...
    t_final, TimeAdapt, TolT, PlotRes, Save, fuel_names, mass_fracL, ...
    inert_comps, mass_fracG, n_saved_solutions, T_0, T_inf, R_0, XRad, ...
    R_end_percent, n_saves, FileName);