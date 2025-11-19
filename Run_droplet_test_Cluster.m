function Run_droplet_test_Cluster(i)
    
    addpath(genpath(pwd)); % Add all subfolders of the current directory to the MATLAB path
    
    Tab      = readtable('Droplet Evaporation.xlsx', 'Sheet', 'Run Sim');
    % nSims    = height(Tab.Numero);
    
    FolderName = num2str(Tab.Numero(i));
    mkdir(FolderName);
    addpath(FolderName);
    
    % ---------------------- MESH AND TIME CONFIGURATION ----------------------
    
    nElems_l          = Tab.nElems_l(i);    % Number of elements in the liquid phase
    nElems_g          = Tab.nElems_g(i);    % Number of elements in the gas phase (first part)
    mu_l              = Tab.mu_l(i);             
    mu_g              = Tab.mu_g(i);    
    p                 = Tab.p(i);           % Degree of the polynomial basis functions
    
    Deltat0           = Tab.Deltat0(i);  % Initial time step
                               % If you get warnings like:
                               %   “Nonlinear solver did not converge at stage 2. Reducing time step”
                               %   or
                               %   "Time error (0.xxxxx) too large. Multiplying time step by a factor 0.xxxxx",
                               % you may reduce Deltat0 to avoid them
                               % However, this is optional — time step control will automatically handle it
    
    t_final           = Tab.t_final(i);     % Final simulation time
                               % Use the actual time you want the simulation to stop at
                               % If you want the simulation to continue until the droplet disappears, 
                               % set a value larger than the expected final time
    
    TimeAdapt         = logical(Tab.TimeAdapt(i));  % Enable adaptive time-stepping (true/false)
    TolT              = Tab.TolT(i);  % Tolerance for temporal error (≥1e-3) (used if TimeAdapt = true)
    
    % ---------------------------- OUTPUT OPTIONS -----------------------------
    
    PlotRes           = logical(Tab.PlotRes(i));   % Plot intermediate results (true/false)
                               % NOTE: if true, Save must be false
    
    Save              = logical(Tab.Save(i));  % Save results to file (true/false)
                               % NOTE: if true, PlotRes must be false
                               
    n_saved_solutions = Tab.n_saved_solutions(i);   % Number of solution snapshots to save throughout simulation
    n_saves           = Tab.n_saves(i);    % Number of intermediate "safety saves" during simulation
    
    % ---------------------- INITIAL AND BOUNDARY CONDITIONS ------------------
    
    R_0               = Tab.R_0(i); % Initial droplet radius [m]
    XRad              = Tab.XRad(i);        % Domain size (in gas phase) as a multiple of the initial radius
    R_end_percent     = Tab.R_end_percent(i);       % Final droplet size as a fraction of initial radius (when Save = true)
    
    T_0               = Tab.T_0(i);        % Initial temperature at the center of the droplet [K]
    T_inf             = Tab.T_inf(i);       % Ambient (far-field) temperature [K]
    
    % --------------------------- SPECIES DEFINITION --------------------------
    
    fuel_names        = strsplit(Tab.fuel_names{i},';');     % Fuel components
    frac_strL         = strsplit(Tab.mass_fracL{i},';');
    mass_fracL        = cellfun(@(x) str2double(strrep(x, ',', '.')), frac_strL);
    mass_fracL        = num2cell(mass_fracL);
    mass_fracL        = mass_fracL(~cellfun(@(x) all(isnan(x)), mass_fracL));   % Mass fraction of each fuel in the liquid phase (same order as above)
    inert_comps       = strsplit(Tab.inert_comps{i},';');    % Inert gas species (DO NOT MODIFY)
    frac_strG         = strsplit(Tab.mass_fracG{i},';'); 
    mass_fracG        = cellfun(@(x) str2double(strrep(x, ',', '.')), frac_strG);
    mass_fracG        = num2cell(mass_fracG);   % Mass fraction of each inert species in the gas phase (same order as above)
    
    % ------------------------ INITIAL SIMULATION CALL ------------------------
    
    save(fullfile(FolderName,'PRUEBA.mat'),'Tab','nElems_l','p','Deltat0','t_final', ...
        'TimeAdapt','TolT','PlotRes','Save','n_saved_solutions','n_saves','R_0','XRad', ...
        'R_end_percent','T_0','T_inf','fuel_names','mass_fracL','inert_comps','mass_fracG')

    % Run the full droplet evaporation simulation from t = 0
    [save_vars] = droplet_test_Cluster2(nElems_l, nElems_g,...
        mu_l, mu_g, p, Deltat0, t_final, TimeAdapt, TolT, ...
        PlotRes, Save, fuel_names, mass_fracL, inert_comps, ...
        mass_fracG, n_saved_solutions, T_0, T_inf, ...
        R_0, XRad, R_end_percent, n_saves, FolderName);
    
end