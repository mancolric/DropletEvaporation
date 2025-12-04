function [save_vars, model_l, model_g] = ...
    droplet_test(nElems_l, nElems_g, hmin_l, hmin_g, p, ...
        Deltat0, t_final, TimeAdapt, TolT, ...
        PlotRes, Save, fuel_names, mass_fracL, ...
        inert_comps, mass_fracG, n_saved_solutions, T_0, T_inf, ...
        x_lg, XRad, R_end_percent, n_saves)

    tStart = tic;

    %----------------------------------------------------------------------
    %DATA:
    
    %Maximum and target nb of iterations in nonlinear solver:
    NLS_MaxIter     = 200;   %for BPR3 method
    NLS_IterTarget  = 120;   %for BPR3 method
    
    %Relaxation time for boundary conditions:
    tau_relax   = 1e-6;
%     tau_relax   = Inf;
    
    % NOTE: a factor controlling the time step is sqrt(NLS_IterTarget/NLS_iters) 
    % where NLS_iters is the mean number of nonlinear solver iterations at 
    % each Runge--Kutta stage. 
    % The value NLS_IterTarget should be adjusted according to the simulation if adaptive
    % time stepping is enabled.
    % This constant controls the update of the number of nonlinear iterations (NLSiters)
    % based on the average value of NLSiters per stage, which is displayed on the results screen,
    % to achieve a stable time step.
    % For example, if it is observed that NLSiters/stage varies approximately between 27 and 30,
    % and many temporary or convergence errors appear in the console, NLS_IterTarget should be a value around 26.
    
    %Characteristic velocities and time:
    NF_vl       = 1e-4;  %Characteristic value for liquid velocity
    NF_vg       = 1e-1;  %Characteristic value for gas velocity
    NF_tau      = 1e-6;  %Characteristic time
    %NOTE: Higher product NF_v*NF_tau: less iterations, less accuracy in the velocity
    
    %Domain limits:
    x_l         = 0.0;
    x_g         = XRad*x_lg;

    %Function to solve hmin*(1+r+r^2+...+r^(N-1))=L:
    function [f,J]=MeshFactor(r, hmin, nElems, L)
        f       = sum(r.^(0:nElems-1)) - L/hmin;
        J       = dot(1:nElems-1, r.^(0:nElems-2));
    end
    %Find mesh factors. Note that the equation above is ill conditioned. It
    %is necessary to start from r_max, which is defined in such a way that 
    %   hmin*r_max^(N-1)=L
    %It is clear that r<r_max, otherwise, hmin(1+r+r^2+...+r^(N-1))>L.
    [r_l, ~, flag]  = NewtonRaphson(...
                        @(r, ComputeJ) MeshFactor(r, hmin_l, nElems_l, x_lg), ...
                        (x_lg/hmin_l)^(1/(nElems_l-1)), 1e-2*x_lg, 0.0, 200);
    if flag<0
        error('Unable to generate mesh for liquid phase')
    end
    [r_g, ~, flag]  = NewtonRaphson(...
                        @(r, ComputeJ) MeshFactor(r, hmin_g, nElems_g, x_g-x_lg), ...
                        ((x_g-x_lg)/hmin_g)^(1/(nElems_g-1)), 1e-2*(x_g-x_lg), 0.0, 200);
    if flag<0
        error('Unable to generate mesh for liquid phase')
    end
    %Generate meshes:
    helems_l        = hmin_l*r_l.^(nElems_l-1:-1:0);
    xmesh_l         = [ 0.0, cumsum(helems_l) ];
    xmesh_l(end)    = x_lg;
    helems_g        = hmin_g*r_g.^(0:nElems_g-1);
    xmesh_g         = x_lg + [ 0.0, cumsum(helems_g) ];
    xmesh_g(end)    = x_g;
    
    %Plot mesh:
    if false
        figure;
        hold on;
        plot(xmesh_l, zeros(size(xmesh_l)), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Liquid'); 
        plot(xmesh_g, zeros(size(xmesh_g)), 'r*', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'Gas'); 
        plot(x_lg, 0, 'go', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Interface');
        xlabel('x [m]');
        title('Point distribution in liquid and gas');
        legend('Location', 'best');
        grid on;
        return
    end

    %Display minimum mesh sizes:
    disp(['hmin_l=', sprintf('%.2E', min(diff(xmesh_l))), ...
            ', hmin_g=', sprintf('%.2E', min(diff(xmesh_g)))])
    
    %Load models:
    model_l        = Liquid_ALE(fuel_names,mass_fracL, inert_comps, mass_fracG);
    model_g        = Gas_ALE(fuel_names, inert_comps, mass_fracG);
    
    %Initial gas and droplet velocity:
    vR_0           = 8e-3;
    w_0            = NaN;   %To be initialized in u0_g
    
    %Initial condition function for the LIQUID phase
    function u = u0_l(x)
              
        T_l         = T_0.*ones(size(x));
        y_l         = mass_fracL;

        rho_l       = calc_rho(model_l,y_l,T_l);
        rhoy_l      = cell(model_l.nSpecies,1);
        for i=1:model_l.nSpecies
            rhoy_l{i} = rho_l.*y_l{i}.*ones(size(x));
        end

        h_l         = calc_h(T_l,y_l,model_l);
        H_l         = rho_l.*h_l;

        v_l         = 0.0*x;

        %Pack all fields into cell arrays
        %u = {ρY1, ..., ρYN, H, v}
        u           = [ rhoy_l; {H_l}; {v_l} ];
        %Each profile is defined as a function of x to match the expected shape 
        %required by the solver: matrices of size (nElems × (2p + 1))

    end

    %Composition at infinity:
    y_gas_inf   = cell(model_g.nInerts,1);
    for i=1:model_g.nInerts
        y_gas_inf{i}   = mass_fracG{i};
    end
    y_fuel_inf  = repmat({0.0}, model_l.nSpecies, 1);
    y_inf       = [y_gas_inf; y_fuel_inf];

    %Initial condition function for the GAS phase
    function u = u0_g_const(x)
        H_g         = cell(1,1);
        rhoy_g      = cell(model_g.nSpecies,1);
        
        T_g         = T_inf + 0*x;
        y_g         = cell(model_g.nSpecies,1);
        for i=1:model_g.nSpecies
            y_g{i}  = y_inf{i} + 0.0*x;
        end

        rho_g       = calc_rho(model_g,y_g,T_g);
        for i=1:model_g.nSpecies
            rhoy_g{i} = rho_g.*y_g{i};
        end
        h_g         = calc_h(T_g,y_g,model_g);
        H_g{1,1}    = h_g.*rho_g;

        v_g         = 0.0*x;
%         v_g         = NF_vg*(x_lg./x).^2;

        %Final vector of initial gas fields
        %u = {ρY1, ..., ρYN, H, v}
        u           = [ rhoy_g; H_g; {v_g} ];
        %Each profile is defined as a function of x to match the expected shape 
        %required by the solver: matrices of size (nElems × (2p + 1))

        %Plot initial condition GAS:
        % figure()
        % plot(MatTranspVec(x),MatTranspVec(y_g{1}), 'DisplayName', 'N2')
        % hold on
        % plot(MatTranspVec(x),MatTranspVec(y_g{2}), 'DisplayName', 'O2')
        % hold on
        % plot(MatTranspVec(x),MatTranspVec(y_g{3}), 'DisplayName', 'CO2')
        % hold on
        % plot(MatTranspVec(x),MatTranspVec(y_g{4}), 'DisplayName', 'H20')
        % hold on
        % plot(MatTranspVec(x),MatTranspVec(y_g{5}), 'DisplayName', 'Heptano')
        % hold on
        % plot(MatTranspVec(x),MatTranspVec(y_g{6}), 'DisplayName', 'Hexadecano')
        % hold on
        % legend('Location', 'best')
        % hold off

        % figure()
        % plot(MatTranspVec(x),MatTranspVec(T_g), 'DisplayName', 'Temperatura Gas')
        % legend('Location', 'best')
    end

    function u = u0_g_slopes(x)
        
        %Left state:
        y_L         = mass_fracL;
        T_L         = T_0;
        rho_L       = calc_rho(model_l,y_L,T_L);
        h_L         = calc_h(T_L,y_L,model_l);
        H_L         = rho_L*h_L;
        
        %Right state:
        T_R         = T_0;
        y_R         = calc_y_gInter(y_L,model_l,model_g,T_R);
        rho_R       = calc_rho(model_g,y_R,T_R);
        [h_R,hi_R]  = calc_h(T_R,y_R,model_g);
        H_R         = rho_R*h_R;
        
        %Mass and thermal diffusion coefficients at the droplet surface (gas side):
        D_T         = calc_D_T(T_R,y_R,model_g);
        D_rho       = calc_D_rho(T_R,y_R,model_g);
        cp_R        = calc_Cp(T_R,y_R,model_g);
        k_R         = rho_R*cp_R*D_T;
        
        %Initial droplet velocity:
        w_0         = -rho_R/(rho_L-rho_R) * vR_0;
        
        %Compute slopes:
        dYR_dx      = zeros(model_g.nSpecies,1);
        sum_Jh      = 0.0;      %Summatory of J_i h_i:
        for ii=1:model_g.nInerts
            fmassR          = 0.0 + 0.0 - rho_R*y_R{ii}*(vR_0-w_0);
            dYR_dx(ii)      = -fmassR/(rho_R*D_rho);
            sum_Jh          = sum_Jh + fmassR*hi_R{ii};
        end
        for ii=model_g.nInerts+1:model_g.nSpecies
            fmassR          = rho_L*y_L{ii-model_g.nInerts}*(0.0-w_0) + ...
                                0.0 - rho_R*y_R{ii}*(vR_0-w_0);
            dYR_dx(ii)      = -fmassR/(rho_R*D_rho);
            sum_Jh          = sum_Jh + fmassR*hi_R{ii};
        end
        fheatR      = H_L*(0.0-w_0) + 0.0 - H_R*(vR_0-w_0);
        dTR_dx      = -(fheatR-sum_Jh)/k_R;
        
        %Create profiles for y_g. Initially, we think on profiles of
        %the form
        %   u = uR + (uinf-uR)*tanh((x-x_lg)/delta)
        %where delta is chosen so that du_dx matches the desired value,
        %i.e.,
        %   (uinf-uR)/delta = du_dx_target
        %If delta>0, i.e., (uinf-uR) and the desired slope have the same
        %sign, the profile goes from uR to uinf. However, if delta<0, the
        %profile will evolve from uR to 2*uR-uinf. Hence, we add a
        %correction term of the form:
        %   (uinf-uR)*(1-sign(delta))*1/2*(tanh((x-x_lg-delta')/delta')+1)
        %with delta' = |delta|.
        y_g         = cell(model_g.nSpecies,1);
        deltav      = zeros(model_g.nSpecies,1);
        for i=1:model_g.nSpecies
            if dYR_dx(i)==0
                delta   = 1.0; %just a random nonzero finite number
                deltap  = 1.0;
            else
                delta   = (y_inf{i}-y_R{i})/dYR_dx(i);
                deltap  = abs(delta);
            end
            y_g{i}      = y_R{i} + (y_inf{i}-y_R{i})*tanh((x-x_lg)/delta) + ...
                            (y_inf{i}-y_R{i})*(1-sign(delta))*...
                            0.5*(tanh((x-x_lg-20*deltap)/deltap)+1);
            deltav(i)   = delta;
        end
        %Exit if any delta is too small:
        if any(abs(deltav)<1e-6)
            display(deltav)
            error(['delta for specie ', num2str(i), ' too small'])
            delta   = sign(delta)*1e-6;
        end
        
        %Ensure the sum of all species mass fractions equals 1 (sanity check)
        aux         = find(deltav<0, 1, 'first');
        if isempty(aux)
            aux     = model_g.nSpecies;
        end
        y_g{aux}        = 0.0*y_g{aux} + 1.0;
        for ii=[1:aux-1,aux+1:model_g.nSpecies]
            y_g{aux}    = y_g{aux} - y_g{ii};
        end
        
        %Create profile for T_g:
        delta       = (T_inf-T_R)/dTR_dx;
        if abs(delta)<1e-6
            error(['delta for temperature too small'])
            delta   = sign(delta)*1e-6;
        end
        deltap      = abs(delta);
        T_g         = T_R + (T_inf-T_R)*tanh((x-x_lg)/delta) + ...
                            (T_inf-T_R)*(1-sign(delta))*...
                            0.5*(tanh((x-x_lg-20*deltap)/deltap)+1) ;
        
        %Compute rhoy and H:
        rho_g       = calc_rho(model_g,y_g,T_g);
        rhoy_g      = cell(model_g.nSpecies,1);
        for i=1:model_g.nSpecies
            rhoy_g{i} = rho_g.*y_g{i};
        end
        h_g         = calc_h(T_g,y_g,model_g);
        H_g         = rho_g.*h_g;

        %Initial GUESS for the velocity (this is not the exact value!):
        v_g         = vR_0 * (x_lg./x).^2;

        %Final vector of initial gas fields:
        %u = {ρY1, ..., ρYN, H, v}
        u           = [ rhoy_g; {H_g}; {v_g} ];
        
%         display(rho_R*CellToVector(y_R))
%         display(H_R)
%         display(dYR_dx)
%         display(dTR_dx)
        
    end

    function u = u0_g_Millan(x)
        
        %Left and right states:
        T_R         = T_0;
        y_L         = mass_fracL;
        y_R         = calc_y_gInter(y_L,model_l,model_g,T_R);
        y_R_fuels   = 0;
        for i=model_g.nInerts+1:model_g.nSpecies
            y_R_fuels = y_R_fuels + y_R{i};
        end
        
        %Initial Condition (A.Millan):
        D_T         = calc_D_T(T_inf,y_inf,model_g);
        T_g         = T_inf+(x_lg./x).*(T_R-T_inf).*erfc((x-x_lg)./(2*sqrt(D_T*1e-5)));
        y_g         = cell(model_g.nSpecies,1);
        for i=1:model_g.nSpecies
            y_g{i}  = y_inf{i}+(x_lg./x).*(y_R{i}-y_inf{i}).*erfc((x-x_lg)./(2*sqrt(D_T*1e-5)));
        end

        %Ensure the sum of all species mass fractions equals 1 (sanity check)
        sum_yG      = zeros(size(x));
        sum_yL      = zeros(size(x));
        for jj=1:model_g.nInerts
            sum_yG  = sum_yG + y_g{jj};
        end
        for jj=model_g.nInerts+1:model_g.nSpecies
            sum_yL  = sum_yL + y_g{jj};
        end
        if max(max(abs((sum_yG+sum_yL)-1)))>1e-6
            disp('ERROR: Mass fractions not equal to 1!')
            return
        end

        rho_g       = calc_rho(model_g,y_g,T_g);
        rhoy_g      = cell(model_g.nSpecies,1);
        for i=1:model_g.nSpecies
            rhoy_g{i} = rho_g.*y_g{i};
        end
        h_g         = calc_h(T_g,y_g,model_g);
        H_g         = h_g.*rho_g;
        v_g         = 0.008286095019622.*((x_lg./x).^2);

        %Final vector of initial gas fields
        %u = {ρY1, ..., ρYN, H, v}
        u           = [ rhoy_g; {H_g}; {v_g} ];
        
    end

    function u = u0_g(x)
%         u       = u0_g_const(x);  
        u       = u0_g_Millan(x);
    end
    
    %Boundary conditions:
    model_l.u1  = @(t) u0_l(x_l);
    model_g.uN  = @(t) u0_g(x_g);

    %Minimum time step:
    Deltat_min  = 1e-12;
                    
    %----------------------------------------------------------------------
    %INITIAL CONDITION (EXACT FOR DIFFERENTIAL VARIABLES, GUESS FOR ALGEBRAIC ONES):

    %Auxiliary variables:
    nDAE_l        = model_l.nDiff+model_l.nAlg;
    nDAE_g        = model_g.nDiff+model_g.nAlg;
    nDiff_lg      = model_g.nDiff;
    
    %Create meshes and finite element spaces:
    mesh_l        = Mesh_Spheric_Create(xmesh_l);
    fes_l         = FES_QX_Create(mesh_l, p);
    fes_lm1       = FES_PX_Create(mesh_l, p-1);
    mesh_g        = Mesh_Spheric_Create(xmesh_g);
    fes_g         = FES_QX_Create(mesh_g, p);
    fes_gm1       = FES_PX_Create(mesh_g, p-1);
    
    %Prolongation matrices for velocity space:
    prolm_l       = PXToQY_Matrix(fes_lm1, fes_l);
    prolm_g       = PXToQY_Matrix(fes_gm1, fes_g);
    
    %Values for PlotFun:
    xiplot        = linspace(-1.0, 1.0, (p+1)^2);
    xplot_l_ini   = PhysicalCoordinates(mesh_l, xiplot);
    xplot_g_ini   = PhysicalCoordinates(mesh_g, xiplot);
    sol_ini_l     = u0_l(xplot_l_ini);
    v_l_ini       = sol_ini_l{model_l.nDiff+1};
    H_l_ini       = sol_ini_l{model_l.nDiff};
    rhoy_l_ini    = sol_ini_l(1:model_l.nSpecies);
    T_l_ini       = calc_T(H_l_ini,rhoy_l_ini,model_l);
    %
    sol_ini_g     = u0_g(xplot_g_ini);
    v_g_ini       = sol_ini_g{model_g.nDiff+1};
    H_g_ini       = sol_ini_g{model_g.nDiff};
    rhoy_g_ini    = sol_ini_g(1:model_g.nSpecies);
    T_g_ini       = calc_T(H_g_ini,rhoy_g_ini,model_g);
    %
    [rho_l_ini,y_l_ini] = calc_rho_y(rhoy_l_ini,model_l);
    [rho_g_ini,y_g_ini] = calc_rho_y(rhoy_g_ini,model_g);
    %
    uplot_l_ini = [ {rho_l_ini}; {T_l_ini}; {v_l_ini}];
    uplot_g_ini = [ {rho_g_ini}; {T_g_ini}; {v_g_ini}];
    
    %Calculate estimated evaporation time (if known, add manually):
    t_evap      = calc_t_evap(T_inf,y_inf,model_l,model_g,mass_fracL,fuel_names,x_lg);
    t_evap      = real(t_evap);

    %Distribution of the saved time instants: 55% are saved in the first 20% of the 
    %simulation and the remaining 45% in the final 80% of the simulation (MODIFIABLE)
    earlySaveFraction = 0.55;
    lateSaveFraction  = 0.45;
    %
    earlyTimePortion  = 0.2;
    lateTimePortion   = 0.8;
    %
    DeltaSaveI        = (t_evap*earlyTimePortion)/(n_saved_solutions*earlySaveFraction);
    DeltaSaveII       = (t_evap*lateTimePortion)/(n_saved_solutions*lateSaveFraction);  

    %Save finite element spaces:
    sol_n.fesl      = fes_l;
    sol_n.fesg      = fes_g;
    sol_n.feslm1    = fes_lm1;
    sol_n.fesgm1    = fes_gm1;
    
    %Set time:
    sol_n.t         = 0.0;
    
    %Initial condition for liquid:
    M_II            = MassMatrix(sol_n.fesl);
    b               = ProjectFun(@u0_l, sol_n.fesl);
    sol_n.ul        = cell(nDAE_l,1);
    sol_n.vl        = cell(model_l.nAlg,1);
    for II=1:model_l.nDiff
        sol_n.ul{II}                = M_II\b((II-1)*sol_n.fesl.nDof+1:II*sol_n.fesl.nDof);
    end
    for II=model_l.nDiff+1:model_l.nDiff+model_l.nAlg
        bProj                       = prolm_l.'*b((II-1)*sol_n.fesl.nDof+1:II*sol_n.fesl.nDof);
        sol_n.vl{II-model_l.nDiff}  = (prolm_l.'*M_II*prolm_l)\bProj;
        sol_n.ul{II}                = prolm_l*sol_n.vl{II-model_l.nDiff};
    end
    
    %Initial condition for gas:
    M_II            = MassMatrix(sol_n.fesg);
    b               = ProjectFun(@u0_g, sol_n.fesg);
    sol_n.ug        = cell(nDAE_g,1);
    sol_n.vg        = cell(model_g.nAlg,1);
    for II=1:model_g.nDiff
        sol_n.ug{II}                = M_II\b((II-1)*sol_n.fesg.nDof+1:II*sol_n.fesg.nDof);
    end
    for II=model_g.nDiff+1:model_g.nDiff+model_g.nAlg
        bProj                       = prolm_g.'*b((II-1)*sol_n.fesg.nDof+1:II*sol_n.fesg.nDof);
        sol_n.vg{II-model_g.nDiff}  = (prolm_g.'*M_II*prolm_g)\bProj;
        sol_n.ug{II}                = prolm_g*sol_n.vg{II-model_g.nDiff};
    end
       
    %Initial values for droplet variables qL, qR, w:
    sol_n.qL    = CellToVector(EvalSolution(sol_n.ul(1:model_l.nDiff), sol_n.fesl, sol_n.fesl.mesh.nElems, 1.0));
    sol_n.qR    = CellToVector(EvalSolution(sol_n.ug(1:model_g.nDiff+model_g.nAlg), sol_n.fesg, 1, -1.0));
    sol_n.w     = 0.0;
    
    %Normalization factors for differential and algebraic variables: 
    NF_l        = max(1e-6, Lqmean(sol_n.ul, sol_n.fesl, 2));
    NF_g        = max(1e-6, Lqmean(sol_n.ug, sol_n.fesg, 2));
    NF_l(end)   = NF_vl;
    NF_g(end)   = NF_vg;
    
    %Set same NF for all densities:
    NF_l(1:model_l.nSpecies)    = max(NF_l(1:model_l.nSpecies));
    NF_g(1:model_g.nSpecies)    = max(NF_g(1:model_g.nSpecies));
    
    %First Save:
    if Save
        folderName = 'results';
        t = sol_n.t;
        if t==0
            sol_n_save    = cell(round((t_evap*earlyTimePortion)/DeltaSaveI+1)+round((t_evap*lateTimePortion*1.3)/DeltaSaveII+1),1);
            sol_n_save{1} = sol_n;

            save(fullfile(folderName,'Model.mat'),'model_l','model_g', 'p', 't_evap')
        end
    end
    
    %Plot during simulation (modificable):
    if PlotRes
        fig1                = figure();
        fig1.Units          = 'normalized';
        fig1.OuterPosition  = [0.5 0 1 1];
    end  

    function PlotFun(sol)
        if PlotRes
            xiplot      = linspace(-1.0, 1.0, (p+1)^2);
            xplot_l     = PhysicalCoordinates(sol.fesl.mesh, xiplot);
            xplot_g     = PhysicalCoordinates(sol.fesg.mesh, xiplot);
            num_sol_l   = EvalSolution(sol.ul, sol.fesl, 1:sol.fesl.mesh.nElems, xiplot);
            num_sol_g   = EvalSolution(sol.ug, sol.fesg, 1:sol.fesg.mesh.nElems, xiplot);
            
            rhoy_l      = cell(model_l.nSpecies,1);
            rhoy_g      = cell(model_g.nSpecies,1);
            for jj = 1:model_l.nSpecies
                rhoy_l{jj} = num_sol_l{jj}; 
            end
            for jj = 1:model_g.nSpecies
                rhoy_g{jj} = num_sol_g{jj}; 
            end
            
            [rho_l,y_l] = calc_rho_y(rhoy_l, model_l);
            [rho_g,~]   = calc_rho_y(rhoy_g, model_g);
            H_l         = num_sol_l{model_l.nDiff};
            H_g         = num_sol_g{model_g.nDiff};
            T_l         = calc_T(H_l,rhoy_l,model_l);
            T_g         = calc_T(H_g,rhoy_g,model_g);
            v_l         = num_sol_l{model_l.nDiff+1};
            v_g         = num_sol_g{model_g.nDiff+1};
            
            uplot_l     = [ {rho_l}; {T_l}; {v_l} ];
            uplot_g     = [ {rho_g}; {T_g}; {v_g} ];
            
            figure(fig1)

            vars        = {'rho_L', 'rho_G', 'T_L', 'T_G', 'v_L', 'v_G', 'w', 'Y_L'};
            for iPlot=1:length(uplot_l)
                subplot(2,length(uplot_l)+1,iPlot)
                hold off
                if iPlot<=2
                    plot(MatTranspVec(xplot_l_ini), MatTranspVec(uplot_l_ini{iPlot}), 'c')
                    hold on
                else
                    plot(MatTranspVec(xplot_l_ini), MatTranspVec(uplot_l_ini{iPlot}), 'w')
                    hold on
                end
                plot(MatTranspVec(xplot_l), MatTranspVec(uplot_l{iPlot}), 'b')
                hold on
                title([vars{2*iPlot-1},', t_n=',num2str(sol.t)])
                subplot(2,length(uplot_l)+1,iPlot+length(uplot_l)+1)
                hold off
                if iPlot<=2
                    plot(MatTranspVec(xplot_g_ini), MatTranspVec(uplot_g_ini{iPlot}), 'm')
                    hold on
                else
                    plot(MatTranspVec(xplot_g_ini), MatTranspVec(uplot_g_ini{iPlot}), 'w')
                    hold on
                end
                plot(MatTranspVec(xplot_g), MatTranspVec(uplot_g{iPlot}), 'r')
                hold on
                title([vars{2*iPlot},', t_n=',num2str(sol.t)])
            end
            
            %Plot mesh velocity:
            [xmesh_l_n, wmesh_l_n, xmesh_g_n, wmesh_g_n] = MeshVelocities(sol);
            subplot(2,length(uplot_l)+1,length(uplot_l)+1)
            hold off
            plot(xmesh_l_n, wmesh_l_n, 'c')
            hold on
            plot(xmesh_g_n, wmesh_g_n, 'm')
            title([vars{end-1},', t_n=',num2str(sol.t)])

            %Plot liquid mass fractions:
            colors = lines(model_l.nSpecies*2);
            subplot(2,length(uplot_l)+1,length(uplot_l)+length(uplot_g)+1+1)
            hold off
            for i=1:model_l.nSpecies
                xlim([0 1.5e-4])
                plot(MatTranspVec(xplot_l_ini), MatTranspVec(y_l_ini{i}), 'Color', colors(i+model_l.nSpecies, :))
                hold on
                plot(MatTranspVec(xplot_l), MatTranspVec(y_l{i}), 'Color', colors(i, :))
                hold on
            end
            title([vars{end},', t_n=',num2str(sol.t)])

            %Plot also boundary conditions:
            num_sol_ql = sol.qL;
            num_sol_qg = sol.qR;

            rhoy_ql      = cell(model_l.nSpecies,1);
            rhoy_qg      = cell(model_g.nSpecies,1);
            for jj = 1:model_l.nSpecies
                rhoy_ql{jj} = num_sol_ql(jj); 
            end
            for jj = 1:model_g.nSpecies
                rhoy_qg{jj} = num_sol_qg(jj); 
            end
            [rho_ql,~] = calc_rho_y(rhoy_ql, model_l);
            [rho_qg,~] = calc_rho_y(rhoy_qg, model_g);
            H_ql       = sol.qL(model_l.nSpecies+1);
            H_qg       = sol.qR(model_g.nSpecies+1);
            T_ql       = calc_T(H_ql,rhoy_ql,model_l);
            T_qg       = calc_T(H_qg,rhoy_qg,model_g);
%             disp(['T_ql=', sprintf('%.8f', T_ql), ', T_qg=', sprintf('%.8f', T_qg)])

            uplot_ql   = [ {rho_ql}; {T_ql} ];
            uplot_qg   = [ {rho_qg}; {T_qg} ];
            for iPlot=1:length(uplot_ql)
                subplot(2,length(uplot_l)+1,iPlot)
                plot(xmesh_l_n(end), uplot_ql{iPlot}, "*b")
                hold on
            end
            for iPlot=1:length(uplot_qg)
                subplot(2,length(uplot_l)+1,iPlot+length(uplot_l)+1)
                plot(xmesh_g_n(1), uplot_qg{iPlot}, "*r")
                hold on
            end

            title_main = '';
            for i = 1:length(fuel_names)
                percentage = mass_fracL{i} * 100; % Convierte a porcentaje
                if i == 1
                    % Para el primer componente
                    title_main = sprintf('%s %.2f%%', fuel_names{i}, percentage);
                else
                    % Para los componentes siguientes, añade un espacio y concatena
                    title_main = sprintf('%s %s %.2f%%', title_main, fuel_names{i}, percentage);
                end
            end
%             sgtitle(title_main);
        end
    end
    PlotFun(sol_n)
            
    %----------------------------------------------------------------------
    %MARCH IN TIME:
    
    %Compute initial residuals DeltaT and DeltaP in coupling conditions for
    %not well-prepared initial conditions:
    y_eq        = [sol_n.qL; sol_n.qR(1:model_g.nDiff)];
    [r,~]       = EquilibriumConditions(model_l, model_g, y_eq, ...
                        0.0, 0.0, NF_l, NF_g, false);
    DeltaT_0    = r(1);
    DeltaP_0    = r(2:1+model_l.nSpecies);

    %The solution at each stage and at t^(n+1) is to be stored in sol_np1:
    sol_np1     = sol_n;
    
    %The unknown is y, which contains the d.o.f. of [rho, H, v, R] at the liquid,
    %[rho, H, v, R] at the gas:
    
    %Number of variables of the blocks in y:
    nNodes_l        = sol_n.fesl.mesh.nElems+1;
    nNodes_g        = sol_n.fesg.mesh.nElems+1;
    N_ul            = model_l.nDiff*sol_n.fesl.nDof;
    N_vl            = model_l.nAlg*sol_n.feslm1.nDof;
    N_ug            = model_g.nDiff*sol_n.fesg.nDof;
    N_vg            = model_g.nAlg*sol_n.fesgm1.nDof;
    N_L             = model_l.nDiff;
    N_R             = nDAE_g;
    N_z             = N_L+N_R+1;
    block_mesh_l    = 1:nNodes_l;
    block_mesh_g    = block_mesh_l(end)+1:block_mesh_l(end)+nNodes_g;
    block_ul        = block_mesh_g(end)+1:block_mesh_g(end)+N_ul;
    block_vl        = block_ul(end)+1:block_ul(end)+N_vl;
    block_ug        = block_vl(end)+1:block_vl(end)+N_ug;
    block_vg        = block_ug(end)+1:block_ug(end)+N_vg;
    block_z         = block_vg(end)+1:block_vg(end)+N_L+N_R+1;
    
    %Prolongation matrix for liquid and gas systems:
    [prolm_l_iv, prolm_l_jv, prolm_l_sv]    = find(prolm_l);
    prolm_ll_iv     = cat(2, 1:N_ul, ...
                            N_ul+prolm_l_iv.');
    prolm_ll_jv     = cat(2, 1:N_ul, ...
                            N_ul+prolm_l_jv.');
    prolm_ll_sv     = cat(2, ones(1,N_ul), prolm_l_sv.');
    prolm_ll        = sparse(prolm_ll_iv, prolm_ll_jv, prolm_ll_sv, ...
                            nDAE_l*sol_n.fesl.nDof, N_ul+N_vl);
    
    [prolm_g_iv, prolm_g_jv, prolm_g_sv]    = find(prolm_g);
    prolm_gg_iv     = cat(2, 1:N_ug, ...
                            N_ug+prolm_g_iv.');
    prolm_gg_jv     = cat(2, 1:N_ug, ...
                            N_ug+prolm_g_jv.');
    prolm_gg_sv     = cat(2, ones(1,N_ug), prolm_g_sv.');
    prolm_gg        = sparse(prolm_gg_iv, prolm_gg_jv, prolm_gg_sv, ...
                            nDAE_g*sol_n.fesg.nDof, N_ug+N_vg);
     
    %Function that evaluates the residual and the Jacobian of the full
    %system of nonlinear equations at each Runge--Kutta stage. 
    %
    %The equations are:
    %   M_f*u_a - b_f_ii - Deltat_n * a_ii * b_f(y_f)       = 0
    %   R_f     - bmesh_f_ii - Deltat_n * a_ii Rdof_f(y_f)  = 0 
    %
    %CAREFUL: This function modifies the variable sol_np1. The values in
    %sol_np1 correspond to the values of y in the last function call.
    function [r,J] = ResidualFun(y, ComputeJ)
        
        %Extract variables: 
        xmesh_l_np1                         = y(block_mesh_l);
        xmesh_g_np1                         = y(block_mesh_g);
        sol_np1.ul(1:model_l.nDiff)         = VectorToCell(y(block_ul), model_l.nDiff);
        sol_np1.vl(1:model_l.nAlg)          = VectorToCell(y(block_vl), model_l.nAlg);
        for II=1:model_l.nAlg
            sol_np1.ul{model_l.nDiff+II}    = prolm_l*sol_np1.vl{II};
        end
        sol_np1.ug(1:model_g.nDiff)         = VectorToCell(y(block_ug), model_g.nDiff);
        sol_np1.vg(1:model_g.nAlg)          = VectorToCell(y(block_vg), model_g.nAlg);
        for II=1:model_g.nAlg
            sol_np1.ug{model_g.nDiff+II}    = prolm_g*sol_np1.vg{II};
        end
        z                   = y(block_z);
        sol_np1.qL          = z(1:N_L);
        sol_np1.qR          = z(N_L+1:N_L+N_R);
        sol_np1.w           = z(end);
        
        %Update matrices for new mesh:
        mesh_l_np1          = Mesh_Spheric_Create(xmesh_l_np1);
        sol_np1.fesl        = FES_QX_Create(mesh_l_np1, p);
        Mm_l_np1            = MassMatrixExpand(MassMatrix(sol_np1.fesl), model_l.nDiff, model_l.nAlg);
        %
        mesh_g_np1          = Mesh_Spheric_Create(xmesh_g_np1);
        sol_np1.fesg        = FES_QX_Create(mesh_g_np1, p);
        Mm_g_np1            = MassMatrixExpand(MassMatrix(sol_np1.fesg), model_g.nDiff, model_g.nAlg);
        
%         PlotFun(sol_np1)
        
        %Update boundary conditions. MATLAB works with copies, not pointers, 
        %so this is not already done in CouplingConditions:
        %left bc for liquid cannot be applied 
        model_l.uN          = @(t) cat(1, VectorToCell(sol_np1.qL, model_l.nDiff));
        model_g.u1          = @(t) cat(1, VectorToCell(sol_np1.qR, nDAE_g));
        
        %------------------------------------------------------------------
        %Update mesh velocity for current solution:
        
        [~, wmesh_l_np1, ~, wmesh_g_np1] = MeshVelocities(sol_np1);
        uw_l                = [ sol_np1.ul; {P1ToQX(wmesh_l_np1, sol_np1.fesl)} ];
        uw_g                = [ sol_np1.ug; {P1ToQX(wmesh_g_np1, sol_np1.fesg)} ];
        
        %------------------------------------------------------------------
        %Residual of mesh motion equation:
        
        kmesh_l_RK(:,ii)    = wmesh_l_np1;
        kmesh_g_RK(:,ii)    = wmesh_g_np1;
        r_mesh_l            = xmesh_l_np1 - bmesh_l_ii - Deltat_n*RKmethod.aI(ii,ii)*kmesh_l_RK(:,ii);
        r_mesh_g            = xmesh_g_np1 - bmesh_g_ii - Deltat_n*RKmethod.aI(ii,ii)*kmesh_g_RK(:,ii);
                        
        %------------------------------------------------------------------
        %IMPOSE DAE FOR LIQUID:
        
        %Compute term due to fluxes and restriction:
        [kDAE_l_RK(:,ii),dfDAE_duw_l,~,dfDAE_dqL,Deltat_CFL_l]   = ...
            FEM_fgQ(model_l, sol_np1.t, uw_l, sol_np1.fesl, ComputeJ);
        
        %Residuals and Jacobians:
        r_l                             = Mm_l_np1*CellToVector(sol_np1.ul) - bDAE_l_ii - RKmethod.aI(ii,ii)*Deltat_n * kDAE_l_RK(:,ii);
        if ComputeJ
            
            %Derivatives of equations related to liquid:
            aux                         = dfDAE_duw_l.jv<=nDAE_l*sol_np1.fesl.nDof; 
            [M_iv, M_jv, M_sv]          = find(Mm_l_np1);
            J_l.iv                      = cat(1, M_iv, dfDAE_duw_l.iv(aux));
            J_l.jv                      = cat(1, M_jv, dfDAE_duw_l.jv(aux));
            J_l.sv                      = cat(1, M_sv, -RKmethod.aI(2,2)*Deltat_n*dfDAE_duw_l.sv(aux));
            
            %Derivatives w.r.t. qL:
            J_lz.iv                     = dfDAE_dqL.iv;
            J_lz.jv                     = dfDAE_dqL.jv;
            J_lz.sv                     = -RKmethod.aI(2,2)*Deltat_n*dfDAE_dqL.sv;
            
            %Derivatives w.r.t. w:
            aux                         = dfDAE_duw_l.jv>nDAE_l*sol_np1.fesl.nDof;
            dfDAE_dwLeg                 = sparse(dfDAE_duw_l.iv(aux), ...
                                                    dfDAE_duw_l.jv(aux)-nDAE_l*sol_np1.fesl.nDof, ...
                                                    dfDAE_duw_l.sv(aux), ...
                                                    nDAE_l*sol_np1.fesl.nDof, ...
                                                    sol_np1.fesl.nDof);
            dwnodes_dw                  = (xmesh_l_np1-xmesh_l_np1(1))/...
                                            (xmesh_l_np1(end)-xmesh_l_np1(1));
            %df_i/dw = df_i/dwLeg_j * dwLeg_j/dwNodes_k * dwNodes_k/dw: 
            dfDAE_dw                    = dfDAE_dwLeg*P1ToQX(dwnodes_dw, sol_np1.fesl);
            J_lz.iv                     = cat(1, J_lz.iv, (1:nDAE_l*sol_np1.fesl.nDof).' );
            J_lz.jv                     = cat(1, J_lz.jv, repmat(N_L+N_R+1, nDAE_l*sol_np1.fesl.nDof, 1));
            J_lz.sv                     = cat(1, J_lz.sv, -RKmethod.aI(2,2)*Deltat_n*dfDAE_dw);
            
        end
        
        %Apply restriction for algebraic dofs:
        r_l                             = prolm_ll.' * r_l;
        if ComputeJ
            %Compute Jacobians and apply restriction:
            J_l_mat                     = sparse(J_l.iv, J_l.jv, J_l.sv, nDAE_l*sol_np1.fesl.nDof, nDAE_l*sol_np1.fesl.nDof);
            J_lz_mat                    = sparse(J_lz.iv, J_lz.jv, J_lz.sv, nDAE_l*sol_np1.fesl.nDof, N_z);
            J_l_mat                     = prolm_ll.' * J_l_mat * prolm_ll;
            J_lz_mat                    = prolm_ll.' * J_lz_mat;
            %Extract values again:
            [J_l.iv, J_l.jv, J_l.sv]    = find(J_l_mat);
            [J_lz.iv, J_lz.jv, J_lz.sv] = find(J_lz_mat);
        end
        
        %Impose null velocity at origin:
        r_l(N_ul+1)                     = y(block_vl(1)) * NF_l(1)*NF_omega_l*Deltat_n; %keep same order of magnitude
        if ComputeJ
            %Delete previous components:
            J_l.sv(J_l.iv==N_ul+1)      = 0.0;
            J_lz.sv(J_lz.iv==N_ul+1)    = 0.0;
            %New components:
            J_l.iv                      = cat(1, J_l.iv, N_ul+1);
            J_l.jv                      = cat(1, J_l.jv, N_ul+1);
            J_l.sv                      = cat(1, J_l.sv, NF_l(1)*NF_omega_l*Deltat_n);
        end
        
%         %Impose known velocity:
%         block_vl                    = model_l.nDiff*sol_np1.fesl.nDof + (1:sol_np1.fesl.nDof).';
%         mass_l                      = MassMatrix(sol_np1.fesl);
%         y_l                         = y(block_l);
%         function v_l=v_l_fun(x)
%             v_l                     = { 0.0*x };
%         end
%         b_vl                        = ProjectFun(@v_l_fun, sol_n.fesl);
%         r_l(block_vl)               = Deltat_n*(mass_l*y_l(block_vl) - b_vl);
%         if ComputeJ
%             %Delete previous components:
%             aux                     = find((J_l.iv>=block_vl(1)) & (J_l.iv<=block_vl(end)));
%             J_l.sv(aux)             = 0.0;
%             %New components:
%             [Mm_iv, Mm_jv, Mm_sv]   = find(mass_l);
%             J_l.iv                  = cat(1, J_l.iv, Mm_iv+block_vl(1)-1);
%             J_l.jv                  = cat(1, J_l.jv, Mm_jv+block_vl(1)-1);
%             J_l.sv                  = cat(1, J_l.sv, Deltat_n*Mm_sv);
%         end 
        
        %------------------------------------------------------------------
        %IMPOSE DAE FOR GAS:
        
        %Compute term due to fluxes and restriction:
        [kDAE_g_RK(:,ii),dfDAE_duw_g,dfDAE_dqR,~,Deltat_CFL_g]   = ...
            FEM_fgQ(model_g, sol_np1.t, uw_g, sol_np1.fesg, ComputeJ);
        
        %Residuals and Jacobians:
        r_g                             = Mm_g_np1*CellToVector(sol_np1.ug) - bDAE_g_ii - RKmethod.aI(ii,ii)*Deltat_n * kDAE_g_RK(:,ii);
        if ComputeJ
            
            %Derivatives of equations related to gas:
            aux                         = dfDAE_duw_g.jv<=nDAE_g*sol_np1.fesg.nDof; %(the contribution of dfDAE_dw*dw_du is ignored for the moment)
            [M_iv, M_jv, M_sv]          = find(Mm_g_np1);
            J_g.iv                      = cat(1, M_iv, dfDAE_duw_g.iv(aux));
            J_g.jv                      = cat(1, M_jv, dfDAE_duw_g.jv(aux));
            J_g.sv                      = cat(1, M_sv, -RKmethod.aI(2,2)*Deltat_n*dfDAE_duw_g.sv(aux));
 
            %Derivatives w.r.t. qR:
            J_gz.iv                     = dfDAE_dqR.iv;
            J_gz.jv                     = N_L + dfDAE_dqR.jv;
            J_gz.sv                     = -RKmethod.aI(2,2)*Deltat_n*dfDAE_dqR.sv;
            
            %Derivatives w.r.t. w:
            aux                         = dfDAE_duw_g.jv>nDAE_g*sol_np1.fesg.nDof;
            dfDAE_dwLeg                 = sparse(dfDAE_duw_g.iv(aux), ...
                                                    dfDAE_duw_g.jv(aux)-nDAE_g*sol_np1.fesg.nDof, ...
                                                    dfDAE_duw_g.sv(aux), ...
                                                    nDAE_g*sol_np1.fesg.nDof, ...
                                                    sol_np1.fesg.nDof);
            dwnodes_dw                  = (xmesh_g_np1-xmesh_g_np1(end))/...
                                            (xmesh_g_np1(1)-xmesh_g_np1(end));
            %df_i/dw = df_i/dwLeg_j * dwLeg_j/dwNodes_k * dwNodes_k/dw: 
            dfDAE_dw                    = dfDAE_dwLeg*P1ToQX(dwnodes_dw, sol_np1.fesg);
            J_gz.iv                     = cat(1, J_gz.iv, (1:nDAE_g*sol_np1.fesg.nDof).' );
            J_gz.jv                     = cat(1, J_gz.jv, repmat(N_L+N_R+1, nDAE_g*sol_np1.fesg.nDof, 1));
            J_gz.sv                     = cat(1, J_gz.sv, -RKmethod.aI(2,2)*Deltat_n*dfDAE_dw);
            
        end
        
        %Apply restriction for algebraic dofs:
        r_g                             = prolm_gg.' * r_g;
        if ComputeJ
            %Compute Jacobians and apply restriction:
            J_g_mat                     = sparse(J_g.iv, J_g.jv, J_g.sv, nDAE_g*sol_np1.fesg.nDof, nDAE_g*sol_np1.fesg.nDof);
            J_gz_mat                    = sparse(J_gz.iv, J_gz.jv, J_gz.sv, nDAE_g*sol_np1.fesg.nDof, N_z);
            J_g_mat                     = prolm_gg.' * J_g_mat * prolm_gg;
            J_gz_mat                    = prolm_gg.' * J_gz_mat;
            %Extract values again:
            [J_g.iv, J_g.jv, J_g.sv]    = find(J_g_mat);
            [J_gz.iv, J_gz.jv, J_gz.sv] = find(J_gz_mat);
        end
        
        %Impose velocity at droplet surface:
        r_g(N_ug+1)                     = (y(block_vg(1))-sol_np1.qR(end)) * NF_g(1)*NF_omega_g*Deltat_n; %keep same order of magnitude
        if ComputeJ
            %Delete previous components:
            J_g.sv(J_g.iv==N_ug+1)      = 0.0;
            J_gz.sv(J_gz.iv==N_ug+1)    = 0.0;
            %New components:
            J_g.iv                      = cat(1, J_g.iv, N_ug+1);
            J_g.jv                      = cat(1, J_g.jv, N_ug+1);
            J_g.sv                      = cat(1, J_g.sv, NF_g(1)*NF_omega_g*Deltat_n);
        end
        
%         %Impose known velocity:
%         block_vg                    = model_g.nDiff*sol_np1.fesg.nDof + (1:sol_np1.fesg.nDof).';
%         mass_g                      = MassMatrix(sol_np1.fesg);
%         y_g                         = y(block_g);
%         function v_g=v_g_fun(x)
%             v_g                     = { sol_np1.qR(end)*(mesh_g_np1.x_faces(1)./x).^2 };
%         end
%         b_vg                        = ProjectFun(@v_g_fun, sol_n.fesg);
%         r_g(block_vg)               = Deltat_n*(mass_g*y_g(block_vg) - b_vg);
%         if ComputeJ
%             %Delete previous components:
%             aux                     = find(J_g.iv>=block_vg(1) & J_g.iv<=block_vg(end));
%             J_g.sv(aux)             = 0.0;
%             %New components:
%             [Mm_iv, Mm_jv, Mm_sv]   = find(mass_g);
%             J_g.iv                  = cat(1, J_g.iv, Mm_iv+block_vg(1)-1);
%             J_g.jv                  = cat(1, J_g.jv, Mm_jv+block_vg(1)-1);
%             J_g.sv                  = cat(1, J_g.sv, Deltat_n*Mm_sv);
%         end 
        
        %------------------------------------------------------------------
        %IMPOSE INTERFACE CONDITIONS:
        
        %Allocate r:
        r_z         = zeros(N_L+N_R+1, 1);

        %Mass and energy fluxes at the left:
        %NOTE: "uwL" is the numerical solution for the liquid, including
        %the mesh velocity
        %"qL" is the Dirichlet condition at the boundary,
        uwL         = EvalSolution(uw_l, sol_np1.fesl, [mesh_l_np1.nElems], 1.0);
        duwL_dx     = EvalSolution_dx(uw_l, sol_np1.fesl, [mesh_l_np1.nElems], 1.0);
        [fL, dfL_duwL, dfL_duwL_dx, dfL_dqL] = ...
            model_l.ftildeN(model_l, sol_np1.t, xmesh_l_np1(end), uwL, duwL_dx, ...
                    (xmesh_l_np1(end)-xmesh_l_np1(end-1))/sol_np1.fesl.p, ComputeJ);
        
        %Mass and energy fluxes at the left:
        %NOTE: "uwL" is the numerical solution for the liquid, including
        %the mesh velocity
        %"qL" is the Dirichlet condition at the boundary,
        uwR         = EvalSolution(uw_g, sol_np1.fesg, [1], -1.0);
        duwR_dx     = EvalSolution_dx(uw_g, sol_np1.fesg, [1], -1.0);
        [fR, dfR_duwR, dfR_duwR_dx, dfR_dqR] = ...
            model_g.ftilde1(model_g, sol_np1.t, xmesh_g_np1(1), uwR, duwR_dx, ...
                    (xmesh_g_np1(2)-xmesh_g_np1(1))/sol_np1.fesg.p, ComputeJ);
                
        %Flux balance for the differential variables:
        r_z(1:nDiff_lg) = [zeros(model_g.nInerts,1); CellToVector(fL)]-...
                            CellToVector(fR); 
        if ComputeJ

            %Allocate:
            J_z.iv      = zeros(0,1);
            J_z.jv      = zeros(0,1);
            J_z.sv      = zeros(0,1);

            %Derivatives w.r.t. qL (Dirichlet conditions at the left of the droplet):
            iv_aux      = zeros(N_L, N_L);
            jv_aux      = zeros(size(iv_aux));
            sv_aux      = zeros(size(iv_aux));
            for II=1:N_L
                for JJ=1:N_L
                    iv_aux(II,JJ)   = model_g.nInerts+II;
                    jv_aux(II,JJ)   = JJ;
                    sv_aux(II,JJ)   = dfL_dqL{II,JJ};
                end
            end
            J_z.iv      = cat(1, J_z.iv, iv_aux(:));
            J_z.jv      = cat(1, J_z.jv, jv_aux(:));
            J_z.sv      = cat(1, J_z.sv, sv_aux(:));

            %Derivatives w.r.t. qR (Dirichlet conditions at the right of the droplet):
            iv_aux      = zeros(model_g.nDiff, N_R);
            jv_aux      = zeros(size(iv_aux));
            sv_aux      = zeros(size(iv_aux));
            for II=1:model_g.nDiff
                for JJ=1:N_R
                    iv_aux(II,JJ)   = II;
                    jv_aux(II,JJ)   = N_L + JJ;
                    sv_aux(II,JJ)   = -dfR_dqR{II,JJ};
                end
            end
            J_z.iv      = cat(1, J_z.iv, iv_aux(:));
            J_z.jv      = cat(1, J_z.jv, jv_aux(:));
            J_z.sv      = cat(1, J_z.sv, sv_aux(:));

            %Derivatives w.r.t. w_droplet are the derivatives w.r.t. w:
            iv_aux      = zeros(N_L, 1);
            jv_aux      = zeros(size(iv_aux));
            sv_aux      = zeros(size(iv_aux));
            for II=1:N_L
                for JJ=1:1
                    iv_aux(II,JJ)   = model_g.nInerts+II;
                    jv_aux(II,JJ)   = N_L + N_R + JJ;
                    sv_aux(II,JJ)   = dfL_duwL{II,end}*1.0 + ...
                                        dfL_duwL_dx{II,end}*(1.0/(xmesh_l_np1(end)-xmesh_l_np1(1)));
                end
            end
            J_z.iv      = cat(1, J_z.iv, iv_aux(:));
            J_z.jv      = cat(1, J_z.jv, jv_aux(:));
            J_z.sv      = cat(1, J_z.sv, sv_aux(:));
            
            %Derivatives w.r.t. w_droplet are the derivatives w.r.t. w:
            iv_aux      = zeros(model_g.nDiff, 1);
            jv_aux      = zeros(size(iv_aux));
            sv_aux      = zeros(size(iv_aux));
            for II=1:model_g.nDiff
                for JJ=1:1
                    iv_aux(II,JJ)   = II;
                    jv_aux(II,JJ)   = N_L + N_R + JJ;
                    sv_aux(II,JJ)   = -dfR_duwR{II,end}*1.0 - ...
                                        dfR_duwR_dx{II,end}*(-1.0/(xmesh_g_np1(end)-xmesh_g_np1(1)));
                end
            end
            J_z.iv      = cat(1, J_z.iv, iv_aux(:));
            J_z.jv      = cat(1, J_z.jv, jv_aux(:));
            J_z.sv      = cat(1, J_z.sv, sv_aux(:));

            %Derivatives w.r.t. uL (numerical solution at the left of the droplet):
            iv_aux      = zeros(N_L, nDAE_l, p+1);
            jv_aux      = zeros(size(iv_aux));
            sv_aux      = zeros(size(iv_aux));
            phim        = sol_np1.fesl.NCompute([1.0]);
            dphim_dxi   = sol_np1.fesl.NCompute([1.0]);
            for II=1:N_L
                for JJ=1:nDAE_l
                    for iDof=1:p+1
                        iv_aux(II,JJ,iDof)  = model_g.nInerts+II;
                        jv_aux(II,JJ,iDof)  = (JJ-1)*sol_np1.fesl.nDof + ...
                                                (mesh_l_np1.nElems-1)*(p+1) + ...
                                                iDof;
                        sv_aux(II,JJ,iDof)  = dfL_duwL{II,JJ}*phim(iDof) + ...
                                                dfL_duwL_dx{II,JJ}*dphim_dxi(iDof)*...
                                                2/(xmesh_l_np1(end)-xmesh_l_np1(end-1));
                    end
                end
            end
            J_zl.iv     = iv_aux(:);
            J_zl.jv     = jv_aux(:);
            J_zl.sv     = sv_aux(:);
            
            %Derivatives w.r.t. uR (numerical solution at the right of the droplet):
            iv_aux      = zeros(nDiff_lg, nDAE_g, p+1);
            jv_aux      = zeros(size(iv_aux));
            sv_aux      = zeros(size(iv_aux));
            phim        = sol_np1.fesg.NCompute([-1.0]);
            dphim_dxi   = sol_np1.fesg.NCompute([-1.0]);
            for II=1:nDiff_lg
                for JJ=1:nDAE_g
                    for iDof=1:p+1
                        iv_aux(II,JJ,iDof)  = II;
                        jv_aux(II,JJ,iDof)  = (JJ-1)*sol_np1.fesg.nDof + ...
                                                iDof;
                        sv_aux(II,JJ,iDof)  = -dfR_duwR{II,JJ}*phim(iDof) - ...
                                                dfR_duwR_dx{II,JJ}*dphim_dxi(iDof)*...
                                                2/(xmesh_g_np1(2)-xmesh_g_np1(1));
                    end
                end
            end
            J_zg.iv     = iv_aux(:);
            J_zg.jv     = jv_aux(:);
            J_zg.sv     = sv_aux(:);
            
        end
        
        %Append equilibrium conditions:
%         DeltaT_np1              = DeltaT_0*exp(-sol_np1.t/tau_relax);
%         DeltaP_np1              = DeltaP_0*exp(-sol_np1.t/tau_relax);
        DeltaT_np1              = DeltaT_0*exp(-(sol_np1.t/tau_relax)^2);
        DeltaP_np1              = DeltaP_0*exp(-(sol_np1.t/tau_relax)^2);
        [rEq, JEq]              = EquilibriumConditions(model_l, model_g, ...
                                    [sol_np1.qL; sol_np1.qR(1:model_g.nDiff)], ...
                                    DeltaT_np1, DeltaP_np1, NF_l, NF_g, true);
        r_z(nDiff_lg+1:end)     = rEq;
        if ComputeJ
            for II=1:model_l.nSpecies+3
                for JJ=1:model_l.nDiff+model_g.nDiff
                    J_z.iv      = cat(1, J_z.iv, nDiff_lg+II);
                    J_z.jv      = cat(1, J_z.jv, JJ);
                    J_z.sv      = cat(1, J_z.sv, JEq(II,JJ));
                end
            end      
        end
        
        %------------------------------------------------------------------
        %OUTPUT:
        
        %Pack variables:
        r       = { r_mesh_l, r_mesh_g, cat(1, r_l, r_g, r_z) };
        if ComputeJ
%             J   = { J_l,    NaN,    J_lz;
%                     NaN,    J_g,    J_gz;
%                     J_zl,   J_zg,   J_z     };
            %Assembly full Jacobian:
            N_l     = N_ul+N_vl;
            N_g     = N_ug+N_vg;
            J.iv    = cat(1, J_l.iv, J_lz.iv, ...
                            J_g.iv+N_l, J_gz.iv+N_l, ...
                            J_zl.iv+N_l+N_g, J_zg.iv+N_l+N_g, J_z.iv+N_l+N_g);
            J.jv    = cat(1, J_l.jv, J_lz.jv+N_l+N_g, ...
                            J_g.jv+N_l, J_gz.jv+N_l+N_g, ...
                            J_zl.jv, J_zg.jv+N_l, J_z.jv+N_l+N_g);
            J.sv    = cat(1, J_l.sv, J_lz.sv, J_g.sv, J_gz.sv, ...
                            J_zl.sv, J_zg.sv, J_z.sv);
            J       = { 1.0, ...
                        1.0, ...
                        sparse(J.iv(:), J.jv(:), J.sv(:), ...
                            N_l+N_g+N_L+N_R+1, N_l+N_g+N_L+N_R+1) };
            
        else 
            J   = NaN;
        end
        
        %Save Deltat for CFL=1 for last stage (first stage is avoided because 
        %velocity at t=0 is unknown):
        if ii==RKmethod.s
            Deltat_CFL_n    = min(Deltat_CFL_l, Deltat_CFL_g);
        end
        
    end

    %Function with the preconditioned residual (necessary for fixed point
    %iteration). Note: We do not solve f(u)=0, but fhat=Sr*f(Sy*yhat)=0, with 
    %Sr and Sy two scaling matrices:
    %The Jacobian is hence Jhat = Sr*df/dy*Sy and we solve
    %g:=Jhat\fhat(yhat)=0
%     fig2    = figure();
    function gscaled = PrecResidualFun(yscaled)
        
        %Compute residual:
        y       = Sy*yscaled;
        [r,~]   = ResidualFun(y, false);   %r={rmesh_l, rmesh_g, [r_l, r_g, r_z]}
        if any(isnan(r{1})) || any(isnan(r{2})) || any(isnan(r{3})) 
            gscaled = NaN;
            return
        end
        
        %Compute preconditioned residual:
        g_mesh_l    = Srinv{1}*r{1};    %Jacobian for this block is the identity
        g_mesh_g    = Srinv{2}*r{2};
        g_lgz       = LUSolve(A_n_fact{3},Srinv{3}*r{3});
        gscaled     = cat(1, g_mesh_l, g_mesh_g, g_lgz);
        
%         figure(fig2)
%         semilogy(abs(g_l), "+-b")
%         hold on
%         block_v     = model_l.nDiff*fes_l.nDof+1:nDAE_l*fes_l.nDof;
%         semilogy(block_v, abs(g_l(block_v)), "x-c")   
%         semilogy(abs(g_g), "+-r")
%         block_v     = model_g.nDiff*fes_g.nDof+1:nDAE_g*fes_g.nDof;
%         semilogy(block_v, abs(g_g(block_v)), "x-m")
%         hold off
        
    end

    %March:
    N_save          = 1;
    FLAG            = 0;
    factor_save     = ((1:n_saves) / n_saves);
    target_times    = factor_save*t_evap;
    times_triggered = false(size(target_times));
    t_prev          = sol_n.t;
    t_save          = DeltaSaveI;
    save_vars       = cell(10,4000);
    Nt              = 0;
    yv_n            = cat(1,    sol_n.fesl.mesh.x_faces, ...
                                sol_n.fesg.mesh.x_faces, ...
                                CellToVector(sol_n.ul(1:model_l.nDiff)), ...
                                CellToVector(sol_n.vl), ...
                                CellToVector(sol_n.ug(1:model_g.nDiff)), ...
                                CellToVector(sol_n.vg), ...
                                sol_n.qL, ...
                                sol_n.qR, ...
                                sol_n.w ) ;
    Deltat_n        = Deltat0;
    Deltat_CFL_n    = NaN;
    while sol_n.t<t_final
        
        %Load RK coefficients:
        if sol_n.t==0.0
            RKmethod    = calcRKmethod_imex('ARS443');
%             RKmethod    = calcRKmethod_imex('BPR3');
%             RKmethod    = calcRKmethod_imex('KC35');
        else
%             RKmethod    = calcRKmethod_imex('BPR3');
            RKmethod    = calcRKmethod_imex('KC35');
        end
        
        %Allocate derivatives:
        kDAE_l_RK       = zeros(nDAE_l*sol_n.fesl.nDof, RKmethod.s);    %For diff variables + alg variables + mesh velocity
        kmesh_l_RK      = zeros(sol_n.fesl.mesh.nElems+1, RKmethod.s);   %For mesh nodes
        kDAE_g_RK       = zeros(nDAE_g*sol_n.fesg.nDof, RKmethod.s);    %For diff variables + alg variables + mesh velocity
        kmesh_g_RK      = zeros(sol_n.fesg.mesh.nElems+1, RKmethod.s);   %For mesh nodes
        
        %Loop until time step has converged:
        RepeatT         = true;
        while RepeatT
            
            %Set time step:
            t_np1       = sol_n.t + Deltat_n;
            if t_np1>=t_final
                t_np1   = t_final;
            end
            Deltat_n    = t_np1-sol_n.t;

            %Characteristic element volumes and maximum radius;
            NF_omega_l      = sum(sol_n.fesl.mesh.V_cells)/sol_n.fesl.mesh.nElems;
            NF_omega_g      = sum(sol_n.fesg.mesh.V_cells)/sol_n.fesg.mesh.nElems;
            R_l             = sol_n.fesl.mesh.x_faces(end);
            R_g             = sol_n.fesg.mesh.x_faces(end);
            
            %Scaling vectors and multiplicity:
            Sy_factors      = [ R_l; R_g; ...
                                NF_l(1:end-1); NF_l(end)*NF_tau/Deltat_n; ...
                                NF_g(1:end-1); NF_g(end)*NF_tau/Deltat_n; ...
                                NF_l(1:end-1); NF_g(1:end); 0.5*(NF_l(end)+NF_g(end)) ];
            Sy_mult         = [ nNodes_l; nNodes_g; ...
                                repmat(sol_n.fesl.nDof, model_l.nDiff, 1); ...
                                repmat(sol_n.feslm1.nDof, model_l.nAlg, 1); ...
                                repmat(sol_n.fesg.nDof, model_g.nDiff, 1); ...
                                repmat(sol_n.fesgm1.nDof, model_g.nAlg, 1); ...
                                repmat(1, N_L, 1); repmat(1, N_R, 1); 1 ];
            Srinv_factors   = 1.0 ./ [ ...
                                NF_omega_l*NF_l(1:end-1); NF_omega_l*NF_l(1)*Deltat_n; ...
                                NF_omega_g*NF_g(1:end-1); NF_omega_g*NF_g(1)*Deltat_n; ...
                                NF_g(end)*NF_g(1:end-1); ...
                                max(DeltaT_0,300); max(DeltaP_0,1e5); NF_l(1); NF_g(1) ];
            Srinv_mult      = [ repmat(sol_n.fesl.nDof, model_l.nDiff, 1);
                                repmat(sol_n.feslm1.nDof, model_l.nAlg, 1); ...
                                repmat(sol_n.fesg.nDof, model_g.nDiff, 1);
                                repmat(sol_n.fesgm1.nDof, model_g.nAlg, 1); ...
                                repmat(1, nDiff_lg, 1); 
                                1; repmat(1, model_l.nSpecies, 1); 1; 1 ]; 
            
            %Scaling matrices:
            Sy              = ScalingMatrix(Sy_factors, Sy_mult);
            Srinv           = { 1.0/R_l * speye(nNodes_l);
                                1.0/R_g * speye(nNodes_g);
                                ScalingMatrix(Srinv_factors, Srinv_mult) }; 
                            
            %Initialize stage variables:
            ii              = 1;
            yv_np1          = yv_n;
            sol_np1         = sol_n;

            %"Solve" first stage. Although the solution is clearly *_ii=*_n, 
            %we need this step to compute the Jacobian:
            bmesh_l_ii      = sol_n.fesl.mesh.x_faces;
            Mm_l_ii         = MassMatrixExpand(MassMatrix(sol_n.fesl), model_l.nDiff, 0);
            My1_l           = cat(1, Mm_l_ii*yv_n(block_ul), zeros(sol_np1.fesl.nDof,1));
            bDAE_l_ii       = My1_l;
            %
            bmesh_g_ii      = sol_n.fesg.mesh.x_faces;
            Mm_g_ii         = MassMatrixExpand(MassMatrix(sol_n.fesg), model_g.nDiff, 0);
            My1_g           = cat(1, Mm_g_ii*yv_n(block_ug), zeros(sol_np1.fesg.nDof,1));
            bDAE_g_ii       = My1_g;
        
            %Save derivatives k(:,1) and Jacobian for first stage:
            [~,A_n]         = ResidualFun(yv_np1, true);
            A_n_fact        = { 1.0; 
                                1.0;
                                LUFactorization(Srinv{3}*A_n{3}*Sy(block_ul(1):end, block_ul(1):end)) };
                                
            %Loop stages:
            NLS_iters       = 0;
            for ii=2:RKmethod.s

                %Update time:
                sol_np1.t           = sol_n.t + Deltat_n*RKmethod.c(ii);
            
                %Auxiliary vectors:
                bmesh_l_ii          = sol_n.fesl.mesh.x_faces + Deltat_n * (kmesh_l_RK(:,1:ii-1)*RKmethod.aI(ii,1:ii-1).');
                bDAE_l_ii           = My1_l + Deltat_n * (kDAE_l_RK(:,1:ii-1)*RKmethod.aI(ii,1:ii-1).');
                bmesh_g_ii          = sol_n.fesg.mesh.x_faces + Deltat_n * (kmesh_g_RK(:,1:ii-1)*RKmethod.aI(ii,1:ii-1).');
                bDAE_g_ii           = My1_g + Deltat_n * (kDAE_g_RK(:,1:ii-1)*RKmethod.aI(ii,1:ii-1).');
                
                %Call nonlinear solver. Derivatives are already saved in kDAE
                %and kmesh. Variables xmesh_ii, wmesh_ii, ..., Mm_ii are also
                %modified according to the values of y. Deltat_CFL_n is
                %also modified when solving the last stage:
                if TimeAdapt && sol_n.t>0
                    TolA            = max(1e-10, 0.01*etaT_nm1);
                else
                    TolA            = 1e-10;
                end
                NDOF                            = nDAE_l*fes_l.nDof + nDAE_g*fes_g.nDof; 
                [yscaled_np1, nIters, NLSFlag]  = Anderson(@PrecResidualFun, Sy\yv_n, ...
                                                        sqrt(NDOF)*TolA, 0.0, NLS_MaxIter, 50);
                yv_np1              = Sy*yscaled_np1;
                if NLSFlag<0
                    break
                end
                NLS_iters           = NLS_iters + nIters/(RKmethod.s-1);
                
%                 %Plot results:
%                 PlotFun(sol_np1)
%                 drawnow()
                
            end
        
            %Save CFL for present time step:
            CFL_n               = Deltat_n/Deltat_CFL_n;
            
            %The considered methods are stiffly accurate, so the
            %solution _ii is the solution at _np1.

            %If NLS has converged, estimate time error and adapt time step
            if NLSFlag>0 && TimeAdapt
                
                %Compute error from embedded RK. The solution at the last
                %stage is given by
                %   M^(n+1) u^(n+1)     = M^1 u^1 + Deltat_n * sum_j b_j k_j
                %The solution of the embedder RK is given by 
                %   M^(n+1) uhat^(n+1)  = M^1 u^1 + Deltat_n * sum_j bhat_j k_j
                %so the error u-uhat is given by
                %   M^(n+1) err^(n+1)   = Deltat_n * sum_j (bhat_j-b) k_j
                b_err           = Deltat_n*(kDAE_l_RK(1:model_l.nDiff*sol_n.fesl.nDof,:)*(RKmethod.bhatI-RKmethod.bI));
                err_np1         = cell(model_l.nDiff, 1);
                Mm_II           = MassMatrix(sol_np1.fesl);
                for II=1:model_l.nDiff
                    err_np1{II} = Mm_II\b_err((II-1)*sol_np1.fesl.nDof+1:II*sol_np1.fesl.nDof);
                end
                etaT_l          = LqNorm(err_np1, NF_l(1:model_l.nDiff), sol_np1.fesl, 2);
                %Idem for gas:
                b_err           = Deltat_n*(kDAE_g_RK(1:model_g.nDiff*sol_np1.fesg.nDof,:)*(RKmethod.bhatI-RKmethod.bI));
                err_np1         = cell(model_g.nDiff, 1);
                Mm_II           = MassMatrix(sol_np1.fesg);
                for II=1:model_g.nDiff
                    err_np1{II} = Mm_II\b_err((II-1)*sol_np1.fesg.nDof+1:II*sol_np1.fesg.nDof);
                end
                etaT_g          = LqNorm(err_np1, NF_g(1:model_g.nDiff), sol_np1.fesg, 2);
                %Total error:
                etaT            = max(etaT_l, etaT_g);
                
                %Apply controler:
                if sol_n.t==0
                    %Go to next time level:
                    Deltat_np1  = Deltat_n;
                    RepeatT     = false;
                elseif etaT<=TolT
                    Deltat_np1  = Deltat_n * min([(0.8*TolT/etaT)^(1.0/RKmethod.order), sqrt(NLS_IterTarget/NLS_iters), 2.0]);
                    RepeatT     = false;
                else
                    rDeltat     = max(0.5, min((0.8*TolT/etaT)^(1.0/RKmethod.order)));
                    Deltat_n    = Deltat_n * rDeltat;
                    RepeatT     = true;
                    disp(['Time error (', num2str(etaT), ') too large. ', ...
                        'Multiplying time step by a factor ', num2str(rDeltat)])
                end
                
            elseif NLSFlag>0 %&&!TimeAdapt
                
                Deltat_np1      = min(1.2*Deltat_n, Deltat0);
                RepeatT         = false;
                etaT            = NaN;
                
            else %decrease time step and try again
                
                disp(['Nonlinear solver did not converge at stage ', num2str(ii),'. Reducing time step'])
                Deltat_n        = Deltat_n*0.5;
                if Deltat_n<Deltat_min
                    sol_np1     = sol_n;
%                     PlotFun(sol_np1)
                    error('Time step too small')
                end
                RepeatT         = true;
            end
            
        end

        %At this point, integration of the time step has converged.
        disp(' ')
        disp([  't=', sprintf('%.2E',sol_n.t),...
                ', tau=', sprintf('%.2E',Deltat_n), ...
                ', CFL=', sprintf('%.2E',CFL_n), ...
                ', errT=', sprintf('%.2E',etaT), ...
                ', NLSiters/stage=', sprintf('%.1f',NLS_iters), ...
                ', tCPU=', sprintf('%.2E', toc(tStart)) ])

        %Evaluate densities and temperature from numerical solution:
        uL               = EvalSolution(sol_np1.ul, sol_n.fesl, [sol_n.fesl.mesh.nElems], [1.0]);
        uR               = EvalSolution(sol_np1.ug, sol_n.fesg, [1], [-1.0]);
        rhoy_ql          = uL(1:model_l.nSpecies);
        rhoy_qg          = uR(1:model_g.nSpecies);
        H_ql             = uL{model_l.nSpecies+1};
        H_qg             = uR{model_g.nSpecies+1};
        [~,y_ql_num]     = calc_rho_y(rhoy_ql, model_l);
        [~,y_qg_num]     = calc_rho_y(rhoy_qg, model_g);
        T_ql_num         = calc_T(H_ql,rhoy_ql,model_l);
        T_qg_num         = calc_T(H_qg,rhoy_qg,model_g);
        %Evaluate densities and temperature from boundary conditions:
        rhoy_ql          = VectorToCell(sol_np1.qL(1:end-1), model_l.nSpecies);
        rhoy_qg          = VectorToCell(sol_np1.qR(1:end-2), model_g.nSpecies);
        H_ql             = sol_np1.qL(model_l.nSpecies+1);
        H_qg             = sol_np1.qR(model_g.nSpecies+1);
        [~,y_ql_bc]      = calc_rho_y(rhoy_ql, model_l);
        [~,y_qg_bc]      = calc_rho_y(rhoy_qg, model_g);
        T_ql_bc          = calc_T(H_ql,rhoy_ql,model_l);
        T_qg_bc          = calc_T(H_qg,rhoy_qg,model_g);
        for II=1:model_g.nInerts
            disp([  'Y_qg=', sprintf('%.4E', 0.0), ...
                    ', DeltaY_ql=', sprintf('%.5E', 0.0), ...
                    ', Y_qg=', sprintf('%.4E', y_qg_bc{II}), ...
                    ', DeltaY_qg=', sprintf('%.5E', abs(y_qg_bc{II}-y_qg_num{II}))])
        end
        for II=model_g.nInerts+1:model_g.nSpecies
            disp([  'Y_qg=', sprintf('%.4E', y_ql_bc{II-model_g.nInerts}), ...
                    ', DeltaY_ql=', sprintf('%.5E', abs(y_ql_bc{II-model_g.nInerts}-y_ql_num{II-model_g.nInerts})), ...
                    ', Y_qg=', sprintf('%.4E', y_qg_bc{II}), ...
                    ', DeltaY_qg=', sprintf('%.5E', abs(y_qg_bc{II}-y_qg_num{II}))])
        end        
        disp([  'T_ql=', sprintf('%.8f', T_ql_bc), ...
                ', DeltaT_ql=', sprintf('%.2E', abs(T_ql_bc-T_ql_num)), ...
                ', T_qg=', sprintf('%.8f', T_qg_bc), ...
                ', DeltaT_qg=', sprintf('%.2E', abs(T_qg_bc-T_qg_num)) ] );
        disp([  'v_l=', sprintf('%.5E', uL{model_l.nDiff+1}), ...
                ', v_g=', sprintf('%.5E', uR{model_g.nDiff+1}), ...
                ', w=', sprintf('%.5E', sol_np1.w) ] );
        
        %Update solution:
        sol_np1.t   = t_np1;    %Correct roundoff errors for last stage
        Nt          = Nt+1;
        Deltat_n    = Deltat_np1;
        sol_n       = sol_np1;
        yv_n        = yv_np1;
        etaT_nm1    = etaT;
        
        %Plot results:
        PlotFun(sol_n)
        drawnow()
        
        %Save results:
        if Save
            num_sol_l  = EvalSolution(sol_n.ul, sol_n.fesl, 1:sol_n.fesl.mesh.nElems, xiplot);
            num_sol_g  = EvalSolution(sol_n.ug, sol_n.fesg, 1:sol_n.fesg.mesh.nElems, xiplot);
    
            rhoy_l     = cell(model_l.nSpecies,1);
            for j = 1:model_l.nSpecies
                rhoy_l{j} = num_sol_l{j}; 
            end
            rhoy_g     = cell(model_g.nSpecies,1);
            for j = 1:model_g.nSpecies
                rhoy_g{j} = num_sol_g{j}; 
            end
            save_rhoy_g = rhoy_g(model_g.nInerts+1:model_g.nSpecies);
            H_l         = num_sol_l{model_l.nDiff};
            T_l         = calc_T(H_l,rhoy_l,model_l);
            save_Tc     = T_l(1,1);
    
            v_l         = num_sol_l{model_l.nDiff+1};
            save_v_l    = v_l(end,end);
            w           = sol_n.w;
    
            mesh_L      = sol_n.fesl.mesh;
            save_d      = mesh_L.x_faces(end,end)*2;
    
            num_sol_ql  = sol_n.qL;
            num_sol_qg  = sol_n.qR;
    
            rhoy_ql     = cell(model_l.nSpecies,1);
            rhoy_qg     = cell(model_g.nSpecies,1);
            for j = 1:model_l.nSpecies
                rhoy_ql{j} = num_sol_ql(j); 
            end
            for j = 1:model_g.nSpecies
                rhoy_qg{j} = num_sol_qg(j); 
            end
            [~,y_ql]         = calc_rho_y(rhoy_ql, model_l);
            [rho_qg,y_qg]    = calc_rho_y(rhoy_qg, model_g);
            save_y_ql        = CellToVector(y_ql);
            save_y_qg        = CellToVector(y_qg);
            H_ql             = sol_n.qL(model_l.nSpecies+1);
            save_Ts          = calc_T(H_ql,rhoy_ql,model_l);
            v_qg             = num_sol_qg(model_g.nSpecies+2);
            save_v_g         = v_qg;
            save_m           = rho_qg(1,1)*(v_qg-w);
    
            save_vars{1,Nt}  = sol_n.t;
            save_vars{2,Nt}  = save_Ts;
            save_vars{3,Nt}  = save_Tc;
            save_vars{4,Nt}  = save_y_ql;
            save_vars{5,Nt}  = save_y_qg;
            save_vars{6,Nt}  = save_m;
            save_vars{7,Nt}  = save_d;
            save_vars{8,Nt}  = save_rhoy_g;
            save_vars{9,Nt}  = save_v_l;
            save_vars{10,Nt} = save_v_g;

            Guide =[
                "Time [s]";
                "Surface temperature [K]";
                "Temperature at the center [K]";
                "Mass fraction in the liquid";
                "Mass fraction in the gas";
                "Evaporated mass flow [kg/s]";
                "Diameter [m]";
                "Density by mass fraction in gas [kg/m3]";
                "Velocity field in the liquid [m/s]";
                "Velocity field in the gas [m/s]"
            ];
    
            t = sol_n.t;
            if t==t_final
                sol_n_save{N_save+1} = sol_n;
                disp('t has arrived to t_final')
                simulationTime       = toc(tStart);
                save(fullfile(folderName, 'Saved_solutionsEnd.mat'),'sol_n_save','simulationTime')
                save(fullfile(folderName, 'Saved_varsEnd.mat'),'save_vars','Guide','simulationTime')
                disp("Simulation Time: " + simulationTime + " [s]")
                return
            elseif sol_n.fesl.mesh.x_faces(end,end)<=x_lg*R_end_percent
                sol_n_save{N_save+1} = sol_n;
                disp(['The surface of the droplet is less than ', sprintf('%.0f', R_end_percent * 100), '% of the surface of the initial one'])
                simulationTime       = toc(tStart);
                save(fullfile(folderName, 'Saved_solutionsEnd.mat'),'sol_n_save','simulationTime')
                save(fullfile(folderName, 'Saved_varsEnd.mat'),'save_vars','Guide','simulationTime')
                disp("Simulation Time: " + simulationTime + " [s]")
                return
            else
                if t<=t_evap*0.2
                    DeltaSave = DeltaSaveI;
                elseif t>t_evap*0.2 && FLAG==0
                    DeltaSave = DeltaSaveII;
                    FLAG      = 1;
                end
                for ii=1:n_saves
                    if ~times_triggered(ii) && (t_prev < target_times(ii)) && (t >= target_times(ii))
                        percentage     = round(factor_save(ii) * 100);
                        Sol_name       = ['Saved_solutions' num2str(percentage) '.mat'];
                        Vars_name      = ['Saved_vars' num2str(percentage) '.mat'];
                        
                        simulationTime = toc(tStart);

                        save(fullfile(folderName, Sol_name), 'sol_n_save','simulationTime');
                        save(fullfile(folderName, Vars_name), 'save_vars','Guide','simulationTime');

                        times_triggered(ii) = true;
                    end
                end
                if (t/t_save)>=1
                    sol_n_save{N_save+1} = sol_n;
                    N_save               = N_save + 1;
                    t_save               = t_save + DeltaSave;
                end
                t_prev = t;
            end
        end
        clearvars sol_np1 RKmethod kDAE_l_RK kmesh_l_RK kDAE_g_RK kmesh_g_RK RepeatT t_np1 Sr_l Sr_g Sy_l Sy_g A_n A_n_fact bmesh_l_ii Mm_l_ii My1_l bDAE_l_ii bmesh_g_ii Mm_g_ii My1_g bDAE_g_ii Deltat_np1 yv_np1 
    end
    
    disp(['Simulation finished, tCPU=', sprintf('%.2E', toc(tStart))])
    
end

%u = [rhobarL, HbarL, rhobarR, HbarR]. We impose
%   equil of temperatures
%   equil of chemical potential
%   density--enthalpy restriction at the liquid
%   density--enthalpy restriction at the gas
function [r,J]  = EquilibriumConditions(model_l, model_g, y_eq, ...
    DeltaT, DeltaP, NF_l, NF_g, ComputeJ)

    rhoy_l        = cell(model_l.nSpecies,1);
    rhoy_g        = cell(model_g.nSpecies,1);
    for II=1:model_l.nSpecies
        rhoy_l{II}  = y_eq(II);
    end
    for II=1:model_g.nSpecies
        rhoy_g{II}  = y_eq(II+model_l.nDiff);
    end

    [rho_l, y_l]  = calc_rho_y(rhoy_l,model_l);
    [rho_g, y_g]  = calc_rho_y(rhoy_g,model_g);

    H_l           = y_eq(model_l.nDiff);
    T_l           = calc_T(H_l,rhoy_l,model_l);
    H_g           = y_eq(model_l.nDiff+model_g.nDiff);
    T_g           = calc_T(H_g,rhoy_g,model_g);

    ec            = calc_ecPvap(y_l, y_g, model_l, model_g, T_l);

    rho_bar_l     = calc_rho(model_l,y_l,T_l);
    rho_bar_g     = calc_rho(model_g,y_g,T_g);

    %Residual and Jacobian:
    r             = [ T_l - T_g - DeltaT;
                      ec - DeltaP;
                      rho_l - rho_bar_l;
                      rho_g - rho_bar_g ];
    if ComputeJ
        J                   = zeros(length(r),length(y_eq));
        NF_l                = NF_l(:);
        NF_g                = NF_g(:);
        NFv                 = cat(1, NF_l(1:model_l.nDiff), NF_g(1:model_g.nDiff));
        for jj=1:length(y_eq)
            delta           = NFv(jj)*1e-6;
            ypert           = y_eq;
            ypert(jj)       = y_eq(jj)-delta;
            [f1,~]          = EquilibriumConditions(model_l, model_g, ypert, DeltaT, DeltaP, NF_l, NF_g, false);
            ypert(jj)       = y_eq(jj)+delta;
            [f2,~]          = EquilibriumConditions(model_l, model_g, ypert, DeltaT, DeltaP, NF_l, NF_g, false);
            J(:,jj)         = (f2-f1)/(2*delta);
        end
    else
        J       = NaN;
    end

end

%Function to extract mesh coordinates and velocities from solution:
function [xmeshl, wmeshl, xmeshg, wmeshg] = MeshVelocities(sol)

    xmeshl      = sol.fesl.mesh.x_faces;
    xmeshg      = sol.fesg.mesh.x_faces;
    wmeshl      = sol.w * (xmeshl-xmeshl(1))/(xmeshl(end)-xmeshl(1));
    wmeshg      = sol.w * (xmeshg-xmeshg(end))/(xmeshg(1)-xmeshg(end));
    
end

%Compute scaling matrix:
function S = ScalingMatrix(alphav, Nv)

    nVars       = length(alphav);
    syv         = zeros(0,1);
    for II=1:nVars
        N       = Nv(II);
        syv     = cat(1, syv, repmat(alphav(II), N, 1));
    end
    Ntotal      = sum(Nv);
    S           = spdiags(syv, 0, Ntotal, Ntotal);
    
end
