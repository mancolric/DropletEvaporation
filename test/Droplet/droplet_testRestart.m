function [save_vars, model_l, model_g] = ...
    droplet_testRestart(~, p, ...
        Deltat0, t_final, TimeAdapt, TolT, ...
        ~, Save, ~, mass_fracL, ...
        ~, mass_fracG, n_saved_solutions, T_0, T_inf, ...
        x_lg, XRad, R_end_percent, n_saves, FileName)

    tStart = tic;

    %----------------------------------------------------------------------
    %DATA:
    
    %Domain limits:
    x_l         = 0.0;
    x_gII       = XRad*x_lg;

    %Load models:
    model       = load('Model.mat');
    model_l     = model.model_l;
    model_g     = model.model_g;
    
    restart         = load(FileName);
    sol_n_save_aux  = restart.sol_n_save;
    sol_n_save_aux  = sol_n_save_aux(~cellfun('isempty', sol_n_save_aux));
    prev_sol        = sol_n_save_aux{end};

    u_l         = prev_sol.ul;
    u_g         = prev_sol.ug;
    fes_l       = prev_sol.fesl;
    fes_g       = prev_sol.fesg;
    a           = prev_sol.w;

    xiplot      = linspace(-1.0, 1.0, p*2+1);

    t_n         = prev_sol.t;
    nElems_l2   = prev_sol.fesl.mesh.nElems;
    nElems_g2   = prev_sol.fesg.mesh.nElems;
    num_sol_l   = EvalSolution(u_l, fes_l, 1:nElems_l2, xiplot);
    num_sol_g   = EvalSolution(u_g, fes_g, 1:nElems_g2, xiplot);
    
    %Restart initial condition:
    function u_restart = u0_l_restart(~)
        rhoy_l      = cell(model_l.nSpecies,1);
        for i = 1:model_l.nSpecies
            rhoy_l{i} = num_sol_l{i}; 
        end
        H_l         = num_sol_l{model_l.nDiff};
        v_l         = num_sol_l{model_l.nDiff+1};
        u_restart = [ rhoy_l; H_l; {v_l} ];
    end

    function u_restart = u0_g_restart(~)
        rhoy_g      = cell(model_g.nSpecies,1);
        for i = 1:model_g.nSpecies
            rhoy_g{i} = num_sol_g{i}; 
        end
        H_g         = num_sol_g{model_g.nDiff};
        v_g         = num_sol_g{model_g.nDiff+1};
        u_restart = [ rhoy_g; H_g; {v_g} ];
    end

    %INITIAL CONDITION:
    %Initial mass fraction in the liquid phase
    function y_l = y0_l(~)
        y_l         = cell(model_l.nSpecies,1);
        for i=1:model_l.nSpecies
            y_l{i} = mass_fracL{i};
        end
    end

    %Initial condition function for the LIQUID phase
    function u = u0_l(x)
        H_l         = cell(1,1);
        rhoy_l      = cell(model_l.nSpecies,1);

        T_l         = T_0.*ones(size(x));

        y_l         = y0_l(x);

        rho_l       = calc_rho(model_l,y_l,T_l);
        for i=1:model_l.nSpecies
            rhoy_l{i} = rho_l.*y_l{i}.*ones(size(x));
        end

        h_l         = calc_h(T_l,y_l,model_l);
        H_l{1}      = h_l.*rho_l;

        v_l         = 0.0*x;
        
        %Pack all fields into cell arrays
        %u = {ρY1, ..., ρYN, H, v}
        u           = [ rhoy_l; H_l; {v_l} ];
        %Each profile is defined as a function of x to match the expected shape 
        %required by the solver: matrices of size (nElems × (2p + 1))
    end

    %Initial condition function for the GAS phase
    function u = u0_g(x)
        H_g         = cell(1,1);
        rhoy_g      = cell(model_g.nSpecies,1);

        T_R         = T_0;
        y_L         = y0_l(x_lg);
        y_R         = calc_y_gInter(y_L,model_l,model_g,T_R);
        y_R_fuels   = 0;
        for i=model_g.nInerts+1:model_g.nSpecies
            y_R_fuels = y_R_fuels + y_R{i};
        end
        y_gas_inf   = cell(model_g.nInerts,1);
        for i=1:model_g.nInerts
            y_gas_inf{i}   = (mass_fracG{i}*0.99999999)';
        end
        y_fuel_inf  = repmat({(1-0.99999999)/model_l.nSpecies}, model_l.nSpecies, 1);
        y_inf       = [y_gas_inf; y_fuel_inf];

        %Initial Condition (D.Betran) [valid ONLY for droplets with up to 2 components]:
        % T_ref       = (2/3)*T_0 + (1/3)*T_inf;
        % y_ref       = cell(model_g.nSpecies,1);
        % for i=1:model_g.nSpecies
        %     y_ref{i} = (2/3)*y_R{i} + (1/3)*y_inf{i};
        % end
        % Dif_g       = calc_D_rho(T_ref,y_ref,model_g);
        % rho_gas     = calc_rho(model_g,y_ref,T_ref);
        % Cp          = calc_Cp(T_ref,y_ref,model_g);
        % Lv_i        = calcula_Lv_fuel(model_l.gota,[],T_ref,[]);
        % Lv          = 0;
        % for i=1:model_l.nSpecies
        %     Lv = Lv + Lv_i(1,i)*y_ref{i+model_g.nInerts};
        % end
        % T_refv      = T_ref(:)';
        % Yi_eval     = zeros(model_g.nInerts, size(T_refv,2)); 
        % for i=1:model_g.nInerts
        % Yi_eval(i,:) = y_ref{i}(:)';
        % end
        % Yf_eval = zeros(model_g.nSpecies-model_g.nInerts, size(T_refv,2));
        % for i=1:model_g.nSpecies-model_g.nInerts
        %     Yf_eval(i,:) = y_ref{i+model_g.nInerts}(:)';
        % end
        % k_g         = MixtureRules('k_gas', T_refv, Yi_eval, Yf_eval, model_g.matrix, model_g.gota, model_g.comp_inerts, model_g.P);
        % 
        % mdot_i      = zeros(model_l.nSpecies,1);
        % mdot        = 4*pi*x_lg*rho_gas*Dif_g.*log(1+y_R_fuels./(1-y_R_fuels));
        % for i=model_g.nInerts+1:model_g.nSpecies
        %     epsilon    = y_R{i}/y_R_fuels;
        %     mdot_i(i-model_g.nInerts)  = mdot*epsilon;
        % end
        % q_sens      = mdot*((Cp*(T_inf-T_R))/(exp((mdot*Cp)/(4*pi()*x_lg*k_g))-1)-Lv);
        % y_g         = cell(model_g.nSpecies,1);
        % for i=1:model_g.nInerts
        %     y_g{i} = (1-(1-exp(-(mdot/(4*pi))./(rho_gas.*Dif_g.*x)))).*frac_masG{i};
        % end
        % for i=model_g.nInerts+1:model_g.nSpecies
        %     y_g{i} = 1-exp(-(mdot_i(i-model_g.nInerts)/(4*pi))./(rho_gas.*Dif_g.*x));
        % end
        % 
        % L_comilla  = Lv+q_sens/mdot;
        % 
        % T_g        = T_R-(L_comilla/Cp)+(T_inf-T_R+(L_comilla/Cp))*exp(-(mdot/(4*pi))*Cp./(k_g.*x));

        %Initial Condition (A.Millan):
        D_T         = calc_D_T(T_inf,y_inf,model_g);
        T_g         = T_inf+(x_lg./x).*(T_R-T_inf).*erfc((x-x_lg)./(2*sqrt(D_T*1e-5)));
        y_g         = cell(model_g.nSpecies,1);
        for i=1:model_g.nSpecies
            y_g{i}      = y_inf{i}+(x_lg./x).*(y_R{i}-y_inf{i}).*erfc((x-x_lg)./(2*sqrt(D_T*1e-5)));
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
        for i=1:model_g.nSpecies
            rhoy_g{i} = rho_g.*y_g{i};
        end
        h_g         = calc_h(T_g,y_g,model_g);
        H_g{1,1}    = h_g.*rho_g;

        v_g         = 0.008286095019622.*((x_lg./x).^2);
        
        %Final vector of initial gas fields
        %u = {ρY1, ..., ρYN, H, v}
        u           = [ rhoy_g; H_g; {v_g} ];
        %Each profile is defined as a function of x to match the expected shape 
        %required by the solver: matrices of size (nElems × (2p + 1))
    end

    %Boundary conditions:
    model_l.u1  = @(t) u0_l(x_l);
    model_g.uN  = @(t) u0_g(x_gII);

    %Minimum time step:
    Deltat_min  = 1e-12;
    
    %Maximum nb of iterations in nonlinear solver:
    NLS_MaxIter     = 200;
    NLS_IterTarget  = 200;
    
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
    NF_tau      = 1e-3;  %Characteristic time
    %NOTE: Higher product NF_v*NF_tau: less iterations, less accuracy in the velocity

    %----------------------------------------------------------------------
    %INITIAL CONDITION (EXACT FOR DIFFERENTIAL VARIABLES, GUESS FOR ALGEBRAIC ONES):

    %Auxiliary variables:
    nDAE_l      = model_l.nDiff+model_l.nAlg;
    nDAE_g      = model_g.nDiff+model_g.nAlg;

    %Calculate estimated evaporation time
    t_evap      = model.t_evap;

    %Distribution of the saved time instants: 55% are saved in the first 20% of the 
    %simulation and the remaining 45% in the final 80% of the simulation (MODIFIABLE)
    earlySaveFraction = 0.55;
    lateSaveFraction  = 0.45;

    earlyTimePortion  = 0.2;
    lateTimePortion   = 0.8;

    DeltaSaveI        = (t_evap*earlyTimePortion)/(n_saved_solutions*earlySaveFraction);
    DeltaSaveII       = (t_evap*lateTimePortion)/(n_saved_solutions*lateSaveFraction);

    %Save finite element spaces:
    sol_n.fesl  = fes_l;
    sol_n.fesg  = fes_g;
    
    %Set time:
    sol_n.t     = t_n;
    
    %Initial condition for liquid:
    M_II        = MassMatrix(sol_n.fesl);
    b           = ProjectFun(@u0_l_restart, sol_n.fesl);
    sol_n.ul    = cell(nDAE_l,1);
    for II=1:nDAE_l
        sol_n.ul{II}    = M_II\b((II-1)*sol_n.fesl.nDof+1:II*sol_n.fesl.nDof);
    end
    
    %Initial condition for gas:
    M_II        = MassMatrix(sol_n.fesg);
    b           = ProjectFun(@u0_g_restart, sol_n.fesg);
    sol_n.ug    = cell(nDAE_g,1);
    for II=1:nDAE_g
        sol_n.ug{II}    = M_II\b((II-1)*sol_n.fesg.nDof+1:II*sol_n.fesg.nDof);
    end
       
    %Initial guess for boundary variables qL, qR:
    sol_n.qL    = CellToVector(EvalSolution(sol_n.ul(1:model_l.nDiff), sol_n.fesl, sol_n.fesl.mesh.nElems, 1.0));
    sol_n.qR    = CellToVector(EvalSolution(sol_n.ug(1:model_g.nDiff+model_g.nAlg), sol_n.fesg, 1, -1.0));
    
    %Initial value for droplet velocity:
    sol_n.w     = a;
    
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
        t = sol_n.t;
        sol_n_save    = cell(round((t_evap*earlyTimePortion)/DeltaSaveI+1)+round((t_evap*lateTimePortion*1.3)/DeltaSaveII+1),1);
        sol_n_save{1} = sol_n;
    end
    
    %----------------------------------------------------------------------
    %MARCH IN TIME:
    
    %The solution at each stage and at t^(n+1) is to be stored in sol_np1:
    sol_np1     = sol_n;
    
    %The unknown is y, which contains the d.o.f. of [rho, H, v] at the liquid,
    %[rho, H, v] at the gas:
    
    %Number of variables of the blocks in y:
    N_l             = nDAE_l*sol_n.fesl.nDof;
    N_g             = nDAE_g*sol_n.fesg.nDof;
    N_L             = model_l.nDiff;
    N_R             = nDAE_g;
    block_l         = 1:N_l;
    block_g         = N_l+1:N_l+N_g;
    block_mesh_l    = block_g(end)+1:block_g(end)+mesh_l.nElems+1;
    block_mesh_g    = block_mesh_l(end)+1:block_mesh_l(end)+mesh_g.nElems+1;
    
    %Function that evaluates the residual and the Jacobian of the full
    %system of nonlinear equations at each Runge--Kutta stage. 
    %
    %The equations are:
    %   M_f*u_a - b_f_ii - Deltat_n * a_ii * b_f(y_f)       = 0
    %   R_f     - bmesh_f_ii - Deltat_n * a_ii Rdof_f(y_f)  = 0 
    %for each phase f.
    %
    %CAREFUL: This function modifies the variable sol_np1. The values in
    %sol_np1 correspond to the values of y in the last function call.
    function [r,J] = ResidualFun(y, ComputeJ)
        
        %Extract variables: 
        y_l                 = y{1};
        sol_np1.ul          = VectorToCell(y_l, nDAE_l);
        %
        y_g                 = y{2};
        sol_np1.ug          = VectorToCell(y_g, nDAE_g);
        xmesh_l_np1         = y{3};
        xmesh_g_np1         = y{4};
        
        %Update matrices for new mesh:
        mesh_l_np1          = Mesh_Spheric_Create(xmesh_l_np1);
        sol_np1.fesl        = FES_QX_Create(mesh_l_np1, p);
        Mm_l_np1            = MassMatrixExpand(MassMatrix(sol_np1.fesl), model_l.nDiff, model_l.nAlg);
        %
        mesh_g_np1          = Mesh_Spheric_Create(xmesh_g_np1);
        sol_np1.fesg        = FES_QX_Create(mesh_g_np1, p);
        Mm_g_np1            = MassMatrixExpand(MassMatrix(sol_np1.fesg), model_g.nDiff, model_g.nAlg);

        %         PlotFun(sol_np1)
        
        %------------------------------------------------------------------
        %Solve coupling conditions:
        
        z0                  = cat(1, sol_np1.qL, sol_np1.qR, sol_np1.w);
        [z, nIters, flag]   = CouplingConditions(z0, sol_np1.t, ...
                                model_l, sol_np1.fesl, sol_np1.ul, NF_l, ...
                                model_g, sol_np1.fesg, sol_np1.ug, NF_g);
        if flag<0
            warning('Unable to solve coupling conditions')
            r               = { NaN, NaN };
            if ComputeJ
                J           = sparse(N_l+N_g, N_l+N_g);
            else
                J           = NaN;
            end
            return
        end
        sol_np1.qL          = z(1:N_L);
        sol_np1.qR          = z(N_L+1:N_L+N_R);
        sol_np1.w           = z(end);
        
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
        [kDAE_l_RK(:,ii),dfDAE_duw_l,~,~,Deltat_CFL_g]   = ...
            FEM_fgQ(model_l, sol_np1.t, uw_l, sol_np1.fesl, ComputeJ);
        
        %Residuals and Jacobians:
        r_l                             = Mm_l_np1*y_l - bDAE_l_ii - RKmethod.aI(ii,ii)*Deltat_n * kDAE_l_RK(:,ii);
        if ComputeJ
            
            %Derivatives of equations related to liquid:
            aux                         = dfDAE_duw_l.jv<=nDAE_l*sol_np1.fesl.nDof; 
            [M_iv, M_jv, M_sv]          = find(Mm_l_np1);
            J_l.iv                      = cat(1, M_iv, dfDAE_duw_l.iv(aux));
            J_l.jv                      = cat(1, M_jv, dfDAE_duw_l.jv(aux));
            J_l.sv                      = cat(1, M_sv, -RKmethod.aI(2,2)*Deltat_n*dfDAE_duw_l.sv(aux));
            
        end
        %         %DEBUG: Force v=0 at the origin:
%         Naux            = sol_np1.fesl.NCompute(-1.0);
%         aux_factor      = 1.0/NF_l(model_l.nDiff+1) * NF_l(1)*mesh_l.V_cells(1)*Deltat_n;
%         r_l(model_l.nDiff*sol_np1.fesl.nDof+1)  = aux_factor*(Naux*sol_np1.ul{model_l.nDiff+1}(1:p+1) - 0.0);
%         if ComputeJ
%             aux             = find(J_l.iv == model_l.nDiff*sol_np1.fesl.nDof+1);
%             J_l.sv(aux)     = 0.0;
%             J_l.iv          = cat(1, J_l.iv, repmat(model_l.nDiff*sol_np1.fesl.nDof+1, p+1,1));
%             J_l.jv          = cat(1, J_l.jv, model_l.nDiff*sol_np1.fesl.nDof + (1:p+1).');
%             J_l.sv          = cat(1, J_l.sv, aux_factor*(Naux(:)));
%         end

        %------------------------------------------------------------------
        %IMPOSE DAE FOR GAS:
        
        %Compute term due to fluxes and restriction:
        [kDAE_g_RK(:,ii),dfDAE_duw_g,~,~,Deltat_CFL_l]   = ...
            FEM_fgQ(model_g, sol_np1.t, uw_g, sol_np1.fesg, ComputeJ);
        
        %Residuals and Jacobians:
        r_g                             = Mm_g_np1*y_g - bDAE_g_ii - RKmethod.aI(ii,ii)*Deltat_n * kDAE_g_RK(:,ii);
        if ComputeJ
            
            %Derivatives of equations related to gas:
            aux                         = dfDAE_duw_g.jv<=nDAE_g*sol_np1.fesg.nDof; %(the contribution of dfDAE_dw*dw_du is ignored for the moment)
            [M_iv, M_jv, M_sv]          = find(Mm_g_np1);
            J_g.iv                      = cat(1, M_iv, dfDAE_duw_g.iv(aux));
            J_g.jv                      = cat(1, M_jv, dfDAE_duw_g.jv(aux));
            J_g.sv                      = cat(1, M_sv, -RKmethod.aI(2,2)*Deltat_n*dfDAE_duw_g.sv(aux));
 
        end

        %------------------------------------------------------------------
        %OUTPUT:
        
        %Pack variables:
        r       = { r_l, r_g, r_mesh_l, r_mesh_g };
        if ComputeJ
            J   = { sparse( J_l.iv, J_l.jv, J_l.sv ), ...
                    sparse( J_g.iv, J_g.jv, J_g.sv ), ...
                    1.0, ...    
                    1.0 };
        else 
            J   = NaN;
        end
        
        %Save Deltat for CFL=1 for last stage (first stage is avoided because 
        %velocity at t=0 is unknown):
        if ii==RKmethod.s
            Deltat_CFL_n    = min(Deltat_CFL_g, Deltat_CFL_l);
        end
        
    end

    %Function with the preconditioned residual (necessary for fixed point
    %iteration). Note: We do not solve f(u)=0, but fhat=Sr*f(Sy*yhat)=0, with 
    %Sr and Sy two scaling matrices:
    %The Jacobian is hence Jhat = Sr*df/dy*Sy and we solve
    %g:=Jhat\fhat(yhat)=0
    function gscaled = PrecResidualFun(yscaled)
        
        %Extract solution at liquid and gas:
        y_l         = Sy_l*yscaled(block_l);
        y_g         = Sy_g*yscaled(block_g);
        y_mesh_l    = Sy_mesh_l*yscaled(block_mesh_l);
        y_mesh_g    = Sy_mesh_g*yscaled(block_mesh_g);
        
        %Compute residual:
        [r,~]   = ResidualFun({y_l, y_g, y_mesh_l, y_mesh_g}, false);   %r={r_l, r_g}
        if any(isnan(r{1})) || any(isnan(r{2}))
            gscaled     = NaN;
            return
        end
        
        %Compute preconditioned residual:
        g_l         = LUSolve(A_n_fact{1},Sr_l*r{1});
        g_g         = LUSolve(A_n_fact{2},Sr_g*r{2});
        g_mesh_l    = Sy_mesh_l\r{3};     %g_mesh = Jhat^{-1} fhat = (Sr I Sy)^{-1} Sr f = Sy^{-1} f
        g_mesh_g    = Sy_mesh_g\r{4};     %g_mesh is dimensionless
        gscaled     = cat(1, g_l, g_g, g_mesh_l, g_mesh_g);
 
    end

    %March:
    N_save          = 1;
    FLAG            = 0;
    factor_save     = ((1:n_saves) / n_saves);
    target_times    = factor_save*t_evap;
    times_triggered = false(size(target_times));
    t_prev          = sol_n.t;
    if t<=t_evap*0.2
        DeltaSave = DeltaSaveI;
    elseif t>t_evap*0.2
        DeltaSave = DeltaSaveII;
    end
    t_save          = t+DeltaSave;
    save_vars       = cell(10,4000);
    Nt              = 0;
    yv_n            = cat(1,    CellToVector(sol_n.ul), ...
                                CellToVector(sol_n.ug), ...
                                sol_n.fesl.mesh.x_faces, ...
                                sol_n.fesg.mesh.x_faces ) ;
    Deltat_n        = Deltat0;
    Deltat_CFL_n    = NaN;
    while sol_n.t<t_final
        
        %Load RK coefficients:
        if Nt==0.0
            RKmethod    = calcRKmethod_imex('ARS443', real(t_evap));
%             RKmethod    = calcRKmethod_imex('EULER');
        else
            RKmethod    = calcRKmethod_imex('BPR3', real(t_evap));
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

            %Scaling matrices for all unkwowns:
            [Sr_l,Sy_l]     = ScalMatrices1(model_l, NF_l, NF_tau, sol_n.fesl, Deltat_n);
            [Sr_g,Sy_g]     = ScalMatrices1(model_g, NF_g, NF_tau, sol_n.fesg, Deltat_n);
            Sy_mesh_l       = sol_n.fesg.mesh.x_faces(end) * speye(mesh_l.nElems+1);
            Sy_mesh_g       = sol_n.fesg.mesh.x_faces(end) * speye(mesh_g.nElems+1);

            %Initialize stage variables:
            ii              = 1;
            yv_np1          = yv_n;
            sol_np1         = sol_n;

            %"Solve" first stage. Although the solution is clearly *_ii=*_n, 
            %we need this step to compute the Jacobian:
            bmesh_l_ii      = sol_n.fesl.mesh.x_faces;
            Mm_l_ii         = MassMatrixExpand(MassMatrix(sol_n.fesl), model_l.nDiff, model_l.nAlg);
            My1_l           = Mm_l_ii*yv_n(block_l);
            bDAE_l_ii       = My1_l;
            %
            bmesh_g_ii      = sol_n.fesg.mesh.x_faces;
            Mm_g_ii         = MassMatrixExpand(MassMatrix(sol_n.fesg), model_g.nDiff, model_g.nAlg);
            My1_g           = Mm_g_ii*yv_n(block_g);
            bDAE_g_ii       = My1_g;
        
            %Save derivatives k(:,1) and Jacobian for first stage:
            [~,A_n]         = ResidualFun({yv_np1(block_l), yv_np1(block_g), ...
                                            yv_np1(block_mesh_l), yv_np1(block_mesh_g)}, ...
                                        true);
            A_n_fact        = { LUFactorization(Sr_l*A_n{1}*Sy_l); 
                                LUFactorization(Sr_g*A_n{2}*Sy_g); 
                                1.0; 
                                1.0 };
                                
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
                %modified according to the values of y.
                if TimeAdapt && sol_n.t>0
                    TolA            = max(1e-8, 0.01*etaT_nm1);
                else
                    TolA            = 1e-8;
                end
                NDOF                            = nDAE_l*fes_l.nDof + nDAE_g*fes_g.nDof; 
                [yscaled_np1, nIters, NLSFlag]  = Anderson(@PrecResidualFun, ...
                                                        cat(1,  Sy_l\yv_np1(block_l), ...
                                                                Sy_g\yv_np1(block_g), ...
                                                                Sy_mesh_l\yv_np1(block_mesh_l), ...
                                                                Sy_mesh_g\yv_np1(block_mesh_g)), ...
                                                        sqrt(NDOF)*TolA, NLS_MaxIter, 50);
                yv_np1                           = cat(1,   Sy_l*yscaled_np1(block_l), ...
                                                            Sy_g*yscaled_np1(block_g), ...
                                                            Sy_mesh_l*yscaled_np1(block_mesh_l), ...
                                                            Sy_mesh_g*yscaled_np1(block_mesh_g) );
                if NLSFlag<0
                    break
                end
                NLS_iters           = NLS_iters + nIters/(RKmethod.s-1);
            end
        
            %Save CFL for present time step:
            CFL_n           = Deltat_n/Deltat_CFL_n;
            
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
%                 etaT            = sqrt(etaT_l^2+etaT_g^2);
                etaT            = max(etaT_l, etaT_g);
                
                %Apply controler:
                if etaT<=TolT
                    % NOTE: sqrt(X/NLS_iters) -> Value depending on the simulation
                    % The value X should be adjusted according to the simulation if adaptive
                    % time stepping is enabled.
                    % This constant controls the update of the number of nonlinear iterations (NLSiters)
                    % based on the average value of NLSiters per stage, which is displayed on the results screen,
                    % to achieve a stable time step.
                    % For example, if it is observed that NLSiters/stage varies approximately between 27 and 30,
                    % and many temporary or convergence errors appear in the console, X should be a value around 26.
                    % It is recommended to adjust this value based on the reported NLSiters/stage results
                    % in the simulation to avoid time steps that are too large or too small.
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
        disp([  't=', sprintf('%.2E',sol_n.t),...
                ', tau=', sprintf('%.2E',Deltat_n), ...
                ', CFL=', sprintf('%.2E',CFL_n), ...
                ', errT=', sprintf('%.2E',etaT), ...
                ', NLSiters/stage=', sprintf('%.1f',NLS_iters), ...
                ', tCPU=', sprintf('%.2E', toc(tStart)) ])
        % disp([  't=', sprintf('%.2E',sol_n.t),...
        %         ', tau=', sprintf('%.2E',Deltat_n), ...
        %         ', NLSiters/stage=', num2str(NLS_iters) ])

        %Update solution:
        Nt          = Nt+1;
        Deltat_n    = Deltat_np1;
        sol_n       = sol_np1;
        yv_n        = yv_np1;
        etaT_nm1    = etaT;
        
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

            folderName = 'results';

            t = sol_n.t;
            if t==t_final
                sol_n_save{N_save+1} = sol_n;
                disp('t has arrived to t_final')
                simulationTime       = toc(tStart);
                save(fullfile(folderName, 'Saved_solutionsEnd.mat'),'sol_n_save','simulationTime')
                save(fullfile(folderName, 'Saved_varsEnd.mat'),'save_vars','simulationTime')
                disp("Simulation Time: " + simulationTime + " [s]")
                return
            elseif sol_n.fesl.mesh.x_faces(end,end)<=x_lg*R_end_percent
                sol_n_save{N_save+1} = sol_n;
                disp(['The surface of the droplet is less than ', sprintf('%.0f', R_end_percent * 100), '% of the surface of the initial one'])
                simulationTime       = toc(tStart);
                save(fullfile(folderName, 'Saved_solutionsEnd.mat'),'sol_n_save','simulationTime')
                save(fullfile(folderName, 'Saved_varsEnd.mat'),'save_vars','simulationTime')
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

                        if ~exist(Sol_name, 'file') && ~exist(Vars_name, 'file')
                            save(fullfile(folderName, Sol_name), 'sol_n_save', 'simulationTime');
                            save(fullfile(folderName, Vars_name), 'save_vars', 'Guide', 'simulationTime');
                    
                            times_triggered(ii) = true;
                            fprintf('Saved: %s and %s\n', Sol_name, Vars_name);
                        end
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
        
end

%u = [rhobarL, HbarL, rhobarR, HbarR]. We impose
%   HbarL           = H_l           (equil of temperatures)
%   HbarR-HbarL     = H_g - H_r     (equil of chemical potential)
%   rhobarL         = rho(HbarL)
%   rhobarR         = rho(HbarR)
function [r,J]  = EquilibriumConditions(model_l, model_g, y_eq, ComputeJ)

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
    r             = [ T_l - T_g;
                      ec;
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
            [f1,~]          = EquilibriumConditions(model_l, model_g, ypert, NF_l, NF_g, false);
            ypert(jj)       = y_eq(jj)+delta;
            [f2,~]          = EquilibriumConditions(model_l, model_g, ypert, NF_l, NF_g, false);
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

%Compute scaling coefficients for variables rhoY_i, H, v:
function [Sr, Sy]   = ScalMatrices1(model, NF, fes, Deltat)

    NF          = NF(:);
    
    %Scale v with NF, multiply algebraic eq. by 1.0/Deltat
    Sy      = ScalingMatrix(    cat(1,  NF(1:end-1), ...
                                        NF(end) ), ...
                                cat(1,  repmat(fes.nDof, model.nDiff+model.nAlg, 1) ) );
%     Sy      = ScalingMatrix(    cat(1,  NF(1:end-1), ...
%                                         NF(end)/Deltat ), ...
%                                 cat(1,  repmat(fes.nDof, model.nDiff+model.nAlg, 1) ) );
    Sr      = ScalingMatrix(    cat(1,  1.0./NF(1:end-1), ...
                                        1.0/NF(1)/Deltat), ...
                                cat(1,  repmat(fes.nDof, model.nDiff+model.nAlg, 1) ) );
                            
end

%Here the vector z=[qL; qR; wdroplet]:
function [z, nIters, flag] = CouplingConditions(z0, t, ...
    model_l, fes_l, u_l, NF_l, ...,
    model_g, fes_g, u_g, NF_g)

    %Extract variables:
    xmesh_l         = fes_l.mesh.x_faces;
    xmesh_g         = fes_g.mesh.x_faces;
    N_L             = model_l.nDiff;
    N_R             = model_g.nDiff+1;
    N_z             = N_L+N_R+1;
    
    %Function with residual and jacobian of restrictions:
    function [r,J]  = ResidualFun(y, ComputeJ)

        %Extract variables:
        qL          = y(1:N_L);
        qR          = y(N_L+1:N_L+N_R);
        w_droplet   = y(end);

        %Allocate r:
        r           = zeros(N_L+N_R+1, 1);

        %Mass and energy fluxes at the left:
        %NOTE: "uwL" is the numerical solution for the liquid, including
        %the mesh velocity
        %"qL" is the Dirichlet condition at the boundary,
        uwL         = [ EvalSolution(u_l, fes_l, [fes_l.mesh.nElems], 1.0); 
                        {w_droplet} ];
        duwL_dx     = [ EvalSolution_dx(u_l, fes_l, [fes_l.mesh.nElems], 1.0);
                        {w_droplet*(1.0/(xmesh_l(end)-xmesh_l(1)))} ];
        model_l.uN  = @(t) cat(1, VectorToCell(qL, N_L));
        [fL, dfL_duwL, dfL_duwL_dx, dfL_dqL] = ...
            model_l.ftildeN(model_l, t, xmesh_l(end), uwL, duwL_dx, ...
                    (xmesh_l(end)-xmesh_l(end-1))/fes_l.p, ComputeJ);

        %Mass and energy fluxes at the right:
        %"uwR" is the numerical solution for the gas, including
        %the velocity
        %"qR" is the Dirichlet condition at the boundary,
        uwR         = [ EvalSolution(u_g, fes_g, 1, -1.0);
                        {w_droplet} ];
        duwR_dx     = [ EvalSolution_dx(u_g, fes_g, 1, -1.0);
                        {w_droplet*(-1.0/(xmesh_g(end)-xmesh_g(1)))} ];
        model_g.u1  = @(t) cat(1, VectorToCell(qR, N_R));
        [fR, dfR_duwR, dfR_duwR_dx, dfR_dqR] = ...
            model_g.ftilde1(model_g, t, xmesh_g(1), uwR, duwR_dx, ...
                    (xmesh_g(2)-xmesh_g(1))/fes_g.p, ComputeJ);

        %Flux balance for the differential variables:
        nDiff_lg        = model_g.nDiff;
        fL              = [zeros(model_g.nInerts,1); CellToVector(fL)];
        r(1:nDiff_lg)   = fL-CellToVector(fR); %nDAE_l=nDAE_g equations
        if ComputeJ

            %Allocate:
            J.iv        = zeros(0,1);
            J.jv        = zeros(0,1);
            J.sv        = zeros(0,1);

            %Derivatives w.r.t. qL (Dirichlet conditions at the left of the droplet):
            iv_aux      = zeros(N_L, N_L);
            jv_aux      = zeros(size(iv_aux));
            sv_aux      = zeros(size(iv_aux));
            for II=1:N_L
                for JJ=1:N_L
                    iv_aux(II,JJ)       = model_g.nInerts+II;
                    jv_aux(II,JJ)       = JJ;
                    sv_aux(II,JJ)       = dfL_dqL{II,JJ};
                end
            end
            J.iv        = cat(1, J.iv, iv_aux(:));
            J.jv        = cat(1, J.jv, jv_aux(:));
            J.sv        = cat(1, J.sv, sv_aux(:));

            %Derivatives w.r.t. qR (Dirichlet conditions at the right of the droplet):
            iv_aux      = zeros(model_g.nDiff, N_R);
            jv_aux      = zeros(size(iv_aux));
            sv_aux      = zeros(size(iv_aux));
            for II=1:model_g.nDiff
                for JJ=1:N_R
                    iv_aux(II,JJ)       = II;
                    jv_aux(II,JJ)       = N_L + JJ;
                    sv_aux(II,JJ)       = -dfR_dqR{II,JJ};
                end
            end
            J.iv        = cat(1, J.iv, iv_aux(:));
            J.jv        = cat(1, J.jv, jv_aux(:));
            J.sv        = cat(1, J.sv, sv_aux(:));

            %Derivatives w.r.t. w_droplet are the derivatives w.r.t. w:
            iv_aux      = zeros(N_L, 1);
            jv_aux      = zeros(size(iv_aux));
            sv_aux      = zeros(size(iv_aux));
            for II=1:N_L
                for JJ=1:1
                    iv_aux(II,JJ)       = model_g.nInerts+II;
                    jv_aux(II,JJ)       = N_L + N_R + JJ;
                    sv_aux(II,JJ)       = dfL_duwL{II,end}*1.0 + ...
                                            dfL_duwL_dx{II,end}*(1.0/(xmesh_l(end)-xmesh_l(1)));
                end
            end
            J.iv        = cat(1, J.iv, iv_aux(:));
            J.jv        = cat(1, J.jv, jv_aux(:));
            J.sv        = cat(1, J.sv, sv_aux(:));
            
            %Derivatives w.r.t. w_droplet are the derivatives w.r.t. w:
            iv_aux      = zeros(model_g.nDiff, 1);
            jv_aux      = zeros(size(iv_aux));
            sv_aux      = zeros(size(iv_aux));
            for II=1:model_g.nDiff
                for JJ=1:1
                    iv_aux(II,JJ)       = II;
                    jv_aux(II,JJ)       = N_L + N_R + JJ;
                    sv_aux(II,JJ)       = -dfR_duwR{II,end}*1.0 - ...
                                            dfR_duwR_dx{II,end}*(-1.0/(xmesh_g(end)-xmesh_g(1)));
                end
            end
            J.iv        = cat(1, J.iv, iv_aux(:));
            J.jv        = cat(1, J.jv, jv_aux(:));
            J.sv        = cat(1, J.sv, sv_aux(:));

        end
 
        %Equilibrium conditions:
        y_eq                    = [qL; qR(1:model_g.nDiff)];
        [rEq, JEq]              = EquilibriumConditions(model_l, model_g, y_eq, NF_l, NF_g, ComputeJ);
        r(nDiff_lg+1:end)       = rEq;
        JEq                     = sparse(JEq);
        [iv,jv,sv]              = find(JEq);
        if ComputeJ
            J.iv                = cat(1, J.iv, nDiff_lg + iv);
            J.jv                = cat(1, J.jv, jv);
            J.sv                = cat(1, J.sv, sv);
        end

        %Construct Jacobian:
        if ComputeJ
            J                   = sparse(J.iv, J.jv, J.sv, N_z, N_z);
        else
            J                   = NaN;
        end
    end
    
    %Scaling matrix:
    Sz      = spdiags( cat(1, NF_l(1:N_L), NF_g(1:N_R), 0.5*(NF_l(N_L+1)+NF_g(N_R))), 0, N_z, N_z);
    
    %Solve coupling conditions:
    [~, J]  = ResidualFun(z0, true);
    J_fact  = LUFactorization(J*Sz);
    function gscaled=PrecondResidual(zscaled)
        z       = Sz*zscaled;
        [r,~]   = ResidualFun(z,false);
        gscaled = LUSolve(J_fact,r);
    end
    [zscaled, nIters, flag]     = Anderson(@PrecondResidual, Sz\z0, 1e-8, 20, 50);
    z                           = Sz*zscaled;
    
end
