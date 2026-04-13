function model=Gas_ALE(fuel_names, comp_inerts, frac_masG)

    %Default values for user-defined parameters:
    model.bool_liq    = false;
    model.fuel        = fuel_names;
    model.comp_inerts = comp_inerts;
    model.Inerts_mw   = [0.0280134 0.031999 0.044009 0.01801528]; %[Kg/mol]
    model.nInerts     = length(model.comp_inerts);
    model.frac_masG   = frac_masG;
    model.species     = [model.comp_inerts , fuel_names];
    model.nSpecies    = length(model.species);
    model.nFuels      = model.nSpecies-model.nInerts;
    model.C           = 1000;
    model.CW          = 50.0;         %Penalty coefficient
    model.R           = 8.3145;       % Universal Gas Cte [J/mol/K]
    % function Q=Qfun(~,x)
    %     Q             = cell(model.nDiff,1);
    %     for II=1:model.nDiff
    %         Q{II}     = 0.0*x;
    %     end
    % end

    function Q=Qfun(model, T, rhoy, x)

        Q             = cell(model.nDiff,1);

        C             = cellfun(@(r, mw) r ./ (1e6*mw), rhoy, num2cell(model.species_mw(:)), 'UniformOutput', false);   %[mol/cm3]
        MatrixExpo    = model.ROrder((1:model.nFuels) + model.nInerts, min(1:model.nSpecies, model.nInerts+1)).';
        CROrder       = cellfun(@(c, exp) abs(c).^exp, repmat(C, 1, model.nFuels), num2cell(MatrixExpo), 'UniformOutput', false);
        
        % for JJ=1:model.nFuels           % JJ = number of reactions
        %     k{JJ}     = model.Arr(JJ+model.nInerts,2)*T.^model.Arr(JJ+model.nInerts ,3).*...
        %                     exp(-model.Arr(JJ+model.nInerts,1)./(model.R.*T));
        %     omega{JJ} = k{JJ}.*prod(cat(3, CROrder{:,JJ}), 3);
        % end

        idx_Arr = (1:model.nFuels) + model.nInerts;
        A1 = model.Arr(idx_Arr, 1);
        A2 = model.Arr(idx_Arr, 2);
        A3 = model.Arr(idx_Arr, 3);
        calc_omega = @(j) 1e-0.*  (A2(j) * T.^A3(j) .* exp(-A1(j) ./ (model.R .* T))) .* prod(cat(3, CROrder{:,j}), 3);
        omega = arrayfun(calc_omega, (1:model.nFuels)', 'UniformOutput', false);
        % for II=2:model.nDiff
        %     for JJ=1:model.nFuels
        %         omega_DiffS{JJ}   = omega{JJ}.*model.DiffStoi(JJ+model.nInerts,II-1); 
        %     end
        %     Q{II} = model.species_mw(II-1) .* sum(cat(3, omega_DiffS{:}), 3);
        % end
        Omega3D = cat(3, omega{:});
        % calc_Q = @(ii) model.species_mw(ii-1) .* sum(Omega3D .* reshape(model.DiffStoi((1:model.nFuels)+model.nInerts, ii-1), 1, 1, []), 3);
        % Q(2:model.nDiff) = arrayfun(calc_Q, (2:model.nDiff)', 'UniformOutput', false);

        calc_Q = @(kk) model.species_mw(kk) .* 1e6 .* sum(Omega3D .* reshape(model.DiffStoi((1:model.nFuels)+model.nInerts, kk), 1, 1, []), 3);
        Q(1:model.nSpecies) = arrayfun(calc_Q, (1:model.nSpecies)', 'UniformOutput', false);
        Q{end}                  = zeros(size(Q{end-1}));
        

        %%%%%%%%%%%%%%%%%%%%%%%% Plots, delete %%%%%%%%%%%%%%%%%%%%%%%%
        figure(3)
        for jj=1:5
            plot(reshape(x'./((250/2)*1e-6),[size(x,2)*size(x,1), 1]), reshape(Q{jj}',[size(Q{jj},2)*size(Q{jj},1), 1]))
            hold on
        end
        hold off
        xlabel("$$r/a_0 \; \left[ - \right]$$", "Interpreter","latex")
        ylabel("$$\dot{m} \, \left[ Kg/(m^3*s) \right]$$", "Interpreter","latex")



        % [~,y_g] = calc_rho_y(rhoy, model);
        % figure(2)
        % 
        % for kk=1:length(y_g)
        %     plot(reshape(x'./((250/2)*1e-6),[size(x,2)*size(x,1), 1]), MatTranspVec(y_g{kk}))
        %     hold on
        % end
        % legend("N2", "O2", "CO2", "H20", "Fuel")
        % hold off
        % xlabel("$$r/a_0 \; \left[ - \right]$$", "Interpreter","latex")
        % ylabel("$$Y_\alpha \, \left[ - \right]$$", "Interpreter","latex")
        % 
        % figure(4)
        % plot(reshape(x'./((250/2)*1e-6),[size(x,2)*size(x,1), 1]), reshape(omega{1}',[size(omega{1},2)*size(omega{1},1), 1]))
        % % xlim([10,11])
        % xlabel("$$r/a_0 \; \left[ - \right]$$", "Interpreter","latex")
        % ylabel("$$\omega \, \left[ mol/(m^3*s) \right]$$", "Interpreter","latex")
        % 
        % figure(5)
        % for kk=1:length(CROrder)
        %     plot(reshape(x'./((250/2)*1e-6),[size(x,2)*size(x,1), 1]), MatTranspVec(CROrder{kk}))
        %     hold on
        % end
        % legend("N2", "O2", "CO2", "H20", "Fuel")
        % hold off
        % xlabel("$$r/a_0 \; \left[ - \right]$$", "Interpreter","latex")
        % ylabel("$$C \, \left[ (mol/cm^3)^{ReacOrder} \right]$$", "Interpreter","latex")
        % 
        % figure(6)
        % for kk=1:length(C)
        %     plot(reshape(x'./((250/2)*1e-6),[size(x,2)*size(x,1), 1]), MatTranspVec(C{kk}))
        %     hold on
        % end
        % legend("N2", "O2", "CO2", "H20", "Fuel")
        % hold off
        % xlabel("$$r/a_0 \; \left[ - \right]$$", "Interpreter","latex")
        % ylabel("$$C \, \left[ mol/cm^3 \right]$$", "Interpreter","latex")

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Q               = repmat({zeros(12, 11)}, 6, 1);
    end
    model.Q           = @Qfun;
    model.u1          = @(t) {    1.0;
                                  1.0; 
                                  1.0;
                                  1.0     };
    model.uN          = @(t) {    1.0;
                                  1.0;
                                  1.0;
                                  NaN     };

    model.gota        = clase_gota(fuel_names);
    model.Lv          = calcula_Lv_fuel(model.gota,[],298.15,[]); % Lv(298.15K)
    model.Hf          = calcula_Hf(model, model.bool_liq);
    model.Arr         = calcula_Arrhenius(model, model.bool_liq);
    model.RStoi       = calcula_R_Stoichiometric(model, model.bool_liq);
    model.ROrder      = calcula_ReactionOrders(model, model.bool_liq);
    model.PStoi       = calcula_P_Stoichiometric(model, model.bool_liq);
    model.DiffStoi    = model.PStoi - model.RStoi;
    model.fuel_mw     = calcula_fuel_mw(model, model.bool_liq); %[Kg/mol]
    model.species_mw  = cat(2,model.Inerts_mw, model.fuel_mw);  %[Kg/mol]
    model.Tmin        = 300; 
    model.Tmax        = 5000;
    model.P           = 101325;
    model.N_polyfit   = 5;
    model.matrix      = Polyfit_properties_pureCompounds(model.gota, model.comp_inerts, model.Tmin, model.Tmax, model.P, model.N_polyfit);
    model.NF          = [];    %Normalization factors
    
    %Mandatory fields:
    model.nDiff       = 1+model.nSpecies;  %Number of differential variables
    model.nAlg        = 1;                 %Number of algebraic variables
    model.nVars       = model.nSpecies+3;  %Total nb of variables, including mesh velocity
    model.fQg         = @fQg;              %Function to compute flux, source terms and restriction
    model.ftilde1     = @ftilde1;          %Function to compute numerical flux at face 1, i.e., impose bondary condition
    model.ftildeN     = @ftildeN;          %Function to compute numerical flux at last face, i.e., impose bondary condition
    model.f_diffusive = @f_diffusive;
    
end
   
%--------------------------------------------------------------------------
%The functions below are only seen through structure model:

%Convection-diffusion flux and source terms:
function [  f, df_du, df_du_dx, ...
            Q, dQ_du, dQ_du_dx, ...
            g, dg_du, dg_du_dx, ...
            lambdav                 ] = ...
    fQg(model, t, x, u, du_dx, ComputeJ)

    %Extract variables:
    [rhoy, ~, rho, ~, y]    = aux_rho(model, u, du_dx);
    H                       = u{model.nDiff};
    v                       = u{end-1};
    w                       = u{end};
    T                       = calc_T(H, rhoy, model);
    
    %Diffusive flux: 
    [f, df_du, df_du_dx, ...
        Dmax]               = f_diffusive(model, u, du_dx, ComputeJ);
    
    %Add convective terms:
    for II=1:model.nDiff
        f{II}               = f{II} + (v-w).*u{II};
        if ComputeJ
            df_du{II,II}    = df_du{II,II} + v-w;
            df_du{II,end-1} = df_du{II,end-1} + u{II};
            df_du{II,end}   = df_du{II,end} - u{II};
        end
    end
    
    %Source:
    Q                       = model.Q(model, T, rhoy, x);
    [dQ_du, dQ_du_dx]       = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, rhoy{1});
    
    %Restriction:
    rho_bar                 = calc_rho(model, y, T);
    g                       = { rho - rho_bar };
    [dg_du, dg_du_dx]       = Cells_Allocate(model.nAlg, model.nVars, ComputeJ, rhoy{1});
    if ComputeJ
        
        %Derivatives w.r.t. rho:
        for II=1:model.nSpecies
            
            %Perturb rhoY:
            delta               = 1e-5;
            rhoY_pert           = rhoy;
            rhoY_pert{II}       = rhoy{II} - delta;
            
            %Evaluate Y_pert, T_pert and rhobar_pert:
            [~, Y_pert]         = calc_rho_y(rhoY_pert, model);
            T_pert              = calc_T(H, rhoY_pert, model);
            rho_bar_pert1       = calc_rho(model, Y_pert, T_pert);
            Q_pert1             = model.Q(model, T, rhoY_pert, x);
            rhoY_pert1          = rhoY_pert;
            
            %Perturb rhoY again:
            rhoY_pert           = rhoy;
            rhoY_pert{II}       = rhoy{II} + delta;
            
            %Evaluate Y_pert, T_pert and rhobar_pert:
            [~, Y_pert]         = calc_rho_y(rhoY_pert, model);
            T_pert              = calc_T(H, rhoY_pert, model);
            rho_bar_pert2       = calc_rho(model, Y_pert, T_pert);
            Q_pert2             = model.Q(model, T, rhoY_pert, x);
            rhoY_pert2          = rhoY_pert;
            
            %Add contribution of rho_bar:
            dg_du{1,II}         = 1.0 - (rho_bar_pert2-rho_bar_pert1)/(2*delta);
            dQ_du(:, II)        = cellfun(@(q2, q1) (q2 - q1) / (2 * delta), Q_pert2, Q_pert1, 'UniformOutput', false);
            
        end
        
        %Derivatives w.r.t. H:
        for II=model.nSpecies+1
            
            %Perturb H:
            delta               = 1e-8;
            H_pert              = H - delta;
            
            %Evaluate rho_bar:
            Y_pert              = y;
            T_pert             = calc_T(H_pert, rhoy, model);
            rho_bar_pert1       = calc_rho(model, Y_pert, T_pert);
            Q_pert1             = model.Q(model, T_pert, rhoy, x);
            
            %Perturb H again:
            H_pert              = H + delta;
            
            %Evaluate rho_bar:
            Y_pert              = y;
            T_pert             = calc_T(H_pert, rhoy, model);
            rho_bar_pert2       = calc_rho(model, Y_pert, T_pert);
            Q_pert2             = model.Q(model, T_pert, rhoy, x);
        
            %Add contribution of rho_bar:
            dg_du{1,II}         = - (rho_bar_pert2-rho_bar_pert1)/(2*delta);
            dQ_du(:, II)        = cellfun(@(q2, q1) (q2 - q1) / (2 * delta), Q_pert2, Q_pert1, 'UniformOutput', false);
            
        end
        
    end
    
    %Maximum Deltat for CFL=1 is of the form:
    %   min ( 1/a0, (h/p)/a1, (h/p)^2/a2 , ...)
    %In particular,
    %   a0  = |dw/dx|,          (mesh distortion)
    %   a1  = |v-w|,            (convection)
    %   a2  = Dmax.             (diffusion)
    %This coefficients are stored in lambdav:
    %Note: a0 could be replaced by |div(v-w)|. However, the latter is
    %slightly more difficult to compute.
    %Vector with maximum characteristic speed and maximum diffusion:
    lambdav                 = { abs(du_dx{end}), abs(v-w), Dmax };
    
end

%Auxiliary function:
function [rhoy, drhoy_dx, rho, sum_drhoy_dx, y] = aux_rho(model, u, du_dx)

    %Extract variables:
    rhoy        = cell(model.nSpecies, 1);
    drhoy_dx    = cell(model.nSpecies, 1);
    for II=1:model.nSpecies
        rhoy{II}     = u{II};
        drhoy_dx{II} = du_dx{II};
    end
   
    rho          = zeros(size(rhoy{1})); % Density of the mixture
    sum_drhoy_dx = zeros(size(rhoy{1}));
    for II=1:model.nSpecies
        rho          = rho + rhoy{II};
        sum_drhoy_dx = sum_drhoy_dx + drhoy_dx{II};
    end
    
    y = cell(model.nSpecies,1);
    for II=1:model.nSpecies
        y{II} = rhoy{II}./rho;
    end

end

%Convective flux:
function [f, df_du] = f_ALE(model, u, ComputeJ)
    
    %Extract velocities:
    v       = u{model.nVars-1};
    w       = u{model.nVars};
    vw      = v-w;
    
    %Compute fluxes:
    f       = cell(model.nDiff, 1);
    df_du   = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, u{1});
    for II=1:model.nDiff
        f{II}                       = u{II}*vw;
        if ComputeJ
            df_du{II,II}            = vw;
            df_du{II,model.nVars-1} = u{II};
            df_du{II,model.nVars}   = -u{II};
        end
    end
    
end

%Diffusive flux:
function [f, df_du, df_du_dx, Dmax] = ...
    f_diffusive(model, u, du_dx, ComputeJ)
    
    %Extract variables:
    [rhoy, drhoy_dx, rho, drho_dx, y] = aux_rho(model, u, du_dx);
    H               = u{model.nDiff};
    dH_dx           = du_dx{model.nDiff};

    %Compute dependent variables:
    T     = calc_T(H, rhoy, model);
    D_T   = calc_D_T(T,y,model);
    D_rho = calc_D_rho(T,y,model);
    h_i   = calc_h_i(T,model);
    
    %Diffusive flux:
    f = cell(model.nDiff,1);
    for II=1:model.nSpecies
        f{II}           = - D_rho.* (drhoy_dx{II} - drho_dx.*y{II});
    end
    f{model.nDiff}      = - D_T.* dH_dx;
    for II=1:model.nSpecies
        f{model.nDiff}  = f{model.nDiff} + D_T.*h_i{II}.*drhoy_dx{II} + f{II}.*h_i{II};
    end
    
    %Allocate derivatives:
    [df_du, df_du_dx]   = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, rhoy{1});
    if ComputeJ
        
        %Mass diffusion flux:
        for II=1:model.nSpecies
            for JJ=1:model.nSpecies
                df_du{II,JJ}                = + D_rho.* drho_dx.* (rho.^-2).* ((II==JJ).*rho - rhoy{II});
            end
        end
        for II=1:model.nSpecies
            for JJ=1:model.nSpecies
                df_du_dx{II,JJ}             = - D_rho.* ((II==JJ) - y{II});
            end
        end

        %Energy flux:
        for II=1:model.nSpecies
            for JJ=1:model.nSpecies
                df_du{model.nSpecies+1,JJ}      = df_du{model.nSpecies+1,JJ} + ...
                                                    df_du{II,JJ}.*h_i{II};
                df_du_dx{model.nSpecies+1,JJ}   = df_du_dx{model.nSpecies+1,JJ} + ...
                                                    D_T.*h_i{II}.*(II==JJ) + ...
                                                    df_du_dx{II,JJ}.*h_i{II};
            end
        end
        df_du_dx{model.nSpecies+1,model.nSpecies+1} = - D_T;

    end

    %Maximum diffusion coefficient:
    Dmax        = max( D_rho, D_T );
    
end

%Penalty flux:
function [f, df_duL, df_duR] = f_penalty(model, uL, uR, hp, ComputeJ)

    %Allocate dummy du_dx:
    du_dx                   = cell(model.nVars, 1);
    for II=1:model.nVars
        du_dx{II}           = NaN*uL{II};
    end
    
    %Diffusion due to penalty:
    function Dv=Dfun(u)
        [rhoy, ~, ~, ~, y]  = aux_rho(model, u, du_dx);
        H                   = u{model.nSpecies+1};
        T                   = calc_T(H, rhoy, model);
        D_T                 = calc_D_T(T,y,model);
        D_rho               = calc_D_rho(T,y,model);
        Dv                  = cat(1, repmat(D_rho, model.nSpecies, 1), D_T);              
    end
    sigmav                  = model.CW*max([Dfun(uL),Dfun(uR)],[],2);
    
    %Set penalty to zero if there is only one specie:
    if model.nSpecies==1
        sigmav(1)           = 0.0;
    end
    
    %Penalty terms:
    f       = cell(model.nDiff, 1);
    df_duL  = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, uL{1});
    df_duR  = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, uR{1});
    for II=1:model.nDiff
        f{II}               = sigmav(II)*(uL{II}-uR{II})./hp;
        if ComputeJ
            df_duL{II,II}   = + sigmav(II)./hp;
            df_duR{II,II}   = - sigmav(II)./hp;
        end
    end
    
end

%df_dq means derivatives of flux w.r.t. the parameters that define the
%boundary condition (rhobar, Hbar, vbar in this case)
function [f, df_du, df_du_dx, df_dq] = ...
    ftilde1(model, t, x, u, du_dx, hp, ComputeJ)
    
    %Left state for convective and penalty terms:
    %   -Differential variables are imposed from boundary conditions.
    %   -Velocity is imposed
    %   -Mesh velocity is extrapolated from solution
    %Diffusive flux is extrapolated.
    u1                  = model.u1(t);                      
    uL                  = cell(model.nVars,1);
    for II=1:model.nDiff
        uL{II}          = u1{II};
    end
    uL{model.nDiff+1}   = u1{model.nDiff+1};
    uL{model.nVars}     = u{model.nVars};
    
    %Right state is the numerical solution:
    uR                  = u;
    duR_dx              = du_dx;    
    
    %Evaluate convective, diffusive and penalty fluxes:
    [fc, dfc_duL]               = f_ALE(model, uL, ComputeJ);
    [fd, dfd_duR, dfd_duR_dx]   = f_diffusive(model, uR, duR_dx, ComputeJ);
    [fp, dfp_duL, dfp_duR]      = f_penalty(model, uL, uR, hp, ComputeJ);
    
    %Compute total flux:
    f           = cell(model.nDiff, 1);
    df_duL      = cell(model.nDiff, model.nVars);
    %df_duL_dx  = 0.0
    df_duR      = cell(model.nDiff, model.nVars);
    df_duR_dx   = cell(model.nDiff, model.nVars);
    for II=1:model.nDiff
        f{II}   = fc{II} + fd{II} + fp{II};
    end
    if ComputeJ
       for II=1:model.nDiff
           for JJ=1:model.nVars
               df_duL{II,JJ}    = dfc_duL{II,JJ} + dfp_duL{II,JJ};
               df_duR{II,JJ}    = dfd_duR{II,JJ} + dfp_duR{II,JJ};
               df_duR_dx{II,JJ} = dfd_duR_dx{II,JJ};
           end
       end
    end
    
    %Derivatives w.r.t. u and q. Recall that 
    %   d/du = d/duL * duL/du + d/duR * duR/du
    %   d/dq = d/duL * duL/dq + d/duR * duR/dq
    df_du           = cell(model.nDiff, model.nVars);
    df_du_dx        = cell(model.nDiff, model.nVars);
    df_dq           = cell(model.nDiff, model.nDiff+model.nAlg);  %Derivatives w.r.t. rhobarY_i, Hbar, vbar
    if ComputeJ
        for II=1:model.nDiff
            for JJ=1:model.nVars
                df_du{II,JJ}    = df_duL{II,JJ}*(JJ==model.nVars) + df_duR{II,JJ}*1.0;
                df_du_dx{II,JJ} = 0.0 + df_duR_dx{II,JJ}*1.0;
            end
            for JJ=1:model.nDiff+model.nAlg %derivative w.r.t. diff vars and velocity
                df_dq{II,JJ}     = df_duL{II,JJ};
            end
        end
    end
    
end

%df_dq means derivatives of flux w.r.t. the parameters that define the
%boundary condition (rhobar, Hbar in this case)
function [f, df_du, df_du_dx, df_dq] = ...
    ftildeN(model, t, x, u, du_dx, hp, ComputeJ)
    
    %Left state is the numerical solution:
    uL              = u;
    duL_dx          = du_dx;    
    
    %Right state:
    %   -Differential variables are imposed from boundary conditions.
    %   -Velocity is extrapolated
    %   -Mesh velocity is extrapolated from solution
    %Diffusive flux is extrapolated:
    uN                  = model.uN(t);                      
    uR                  = cell(model.nVars,1);
    for II=1:model.nDiff
        uR{II}          = uN{II};
    end
    uR{model.nDiff+1}   = u{model.nDiff+1};
    uR{model.nVars}     = u{model.nVars};
    
    %Evaluate convective, diffusive and penalty fluxes:
    [fc, dfc_duL]               = f_ALE(model, uL, ComputeJ);
    [fd, dfd_duL, dfd_duL_dx]   = f_diffusive(model, uL, duL_dx, ComputeJ);
    [fp, dfp_duL, dfp_duR]      = f_penalty(model, uL, uR, hp, ComputeJ);
    
    %Compute total flux:
    f           = cell(model.nDiff, 1);
    df_duL      = cell(model.nDiff, model.nVars);
    df_duL_dx   = cell(model.nDiff, model.nVars);
    df_duR      = cell(model.nDiff, model.nVars);
    %df_duR_dx  = 0.0
    for II=1:model.nDiff
        f{II}   = fc{II} + fd{II} + fp{II};
    end
    if ComputeJ
       for II=1:model.nDiff
           for JJ=1:model.nVars
               df_duL{II,JJ}    = dfc_duL{II,JJ} + dfd_duL{II,JJ} +  dfp_duL{II,JJ};
               df_duR{II,JJ}    = dfp_duR{II,JJ};
               df_duL_dx{II,JJ} = dfd_duL_dx{II,JJ};
           end
       end
    end
    
    %Derivatives w.r.t. u and q. Recall that 
    %   d/du = d/duL * duL/du + d/duR * duR/du
    %   d/dq = d/duL * duL/dq + d/duR * duR/dq
    df_du           = cell(model.nDiff, model.nVars);
    df_du_dx        = cell(model.nDiff, model.nVars);
    df_dq           = cell(model.nDiff, model.nDiff);  %Derivatives w.r.t. rhobarY_i, Hbar
    if ComputeJ
        for II=1:model.nDiff
            for JJ=1:model.nVars
                df_du{II,JJ}    = df_duL{II,JJ}*1.0 + df_duR{II,JJ}*(JJ==model.nVars-1) + df_duR{II,JJ}*(JJ==model.nVars);
                df_du_dx{II,JJ} = df_duL_dx{II,JJ}*1.0 + 0.0;
            end
            for JJ=1:model.nDiff %derivative w.r.t. diff vars
                df_dq{II,JJ}    = df_duR{II,JJ};
            end
        end
    end
    
end
