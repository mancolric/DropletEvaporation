function model=Gas_ALE(fuel_names, comp_inerts, frac_masG)

    %Default values for user-defined parameters:
    model.bool_liq    = false;
    model.fuel        = fuel_names;
    model.comp_inerts = comp_inerts;
    model.Inerts_mw   = [0.0280134 0.031999 0.044009 0.01801528];
    model.nInerts     = length(model.comp_inerts);
    model.frac_masG   = frac_masG;
    model.species     = [model.comp_inerts , fuel_names];
    model.nSpecies    = length(model.species);
    model.C           = 1000;
    model.CW          = 50.0;         %Penalty coefficient
    function Q=Qfun(~,x)
        Q             = cell(model.nDiff,1);
        for II=1:model.nDiff
            Q{II}     = 0.0*x;
        end
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
    model.Tmin        = 300; 
    model.Tmax        = 1700;
    model.P           = 101325;
    model.N_polyfit   = 5;
    model.matrix      = Polyfit_properties_pureCompounds(model.gota, model.comp_inerts, model.Tmin, model.Tmax, model.P, model.N_polyfit);
    model.tau_g       = Inf;   %Characteristic time for stabilization restriction; tau_g=Inf disables the latter
    model.NF          = [];    %Normalization factors
    
    %Mandatory fields:
    model.nDiff       = 1+model.nSpecies;  %Number of differential variables
    model.nAlg        = 1;                 %Number of algebraic variables
    model.nVars       = model.nSpecies+3;  %Total nb of variables, including mesh velocity
    model.fQg         = @fQg;              %Function to compute flux, source terms and restriction
    model.ftilde      = @ftilde;           %Function to compute numerical flux at the internal faces
    model.ftilde1     = @ftilde1;          %Function to compute numerical flux at face 1, i.e., impose bondary condition
    model.ftildeN     = @ftildeN_Dirichlet;%Function to compute numerical flux at last face, i.e., impose bondary condition
%     model.ftildeN     = @ftildeN_NoFlux;    %Function to compute numerical flux at last face, i.e., impose bondary condition
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
    
    
    rho_bar = calc_rho(model, y, T);

    %Source:
    Q                       = model.Q(t,x);
    [dQ_du, dQ_du_dx]       = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, rhoy{1});
    
    %Restriction:
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
            
            %Perturb rhoY again:
            rhoY_pert           = rhoy;
            rhoY_pert{II}       = rhoy{II} + delta;
            
            %Evaluate Y_pert, T_pert and rhobar_pert:
            [~, Y_pert]         = calc_rho_y(rhoY_pert, model);
            T_pert              = calc_T(H, rhoY_pert, model);
            rho_bar_pert2       = calc_rho(model, Y_pert, T_pert);
            
            %Add contribution of rho_bar:
            dg_du{1,II}         = 1.0 - (rho_bar_pert2-rho_bar_pert1)/(2*delta);
            
        end
        
        %Derivatives w.r.t. H:
        for II=model.nSpecies+1
            
            %Perturb H:
            delta               = 1e-4;
            H_pert              = H - delta;
            
            %Evaluate rho_bar:
            Y_pert              = y;
            T_pert              = calc_T(H_pert, rhoy, model);
            rho_bar_pert1       = calc_rho(model, Y_pert, T_pert);
            
            %Perturb H again:
            H_pert              = H + delta;
            
            %Evaluate rho_bar:
            Y_pert              = y;
            T_pert              = calc_T(H_pert, rhoy, model);
            rho_bar_pert2       = calc_rho(model, Y_pert, T_pert);
        
            %Add contribution of rho_bar:
            dg_du{1,II}         = - (rho_bar_pert2-rho_bar_pert1)/(2*delta);
            
        end
        
    end
    
    %Add stabilization for density:
    for II=1:model.nSpecies
        Q{II}               = Q{II} - 1/model.tau_g * y{II}.*g{1};
    end
    if ComputeJ
        for II=1:model.nSpecies
            for JJ=1:model.nDiff
                dQ_du{II,JJ}    = dQ_du{II,JJ} - 1/model.tau_g * (...
                                    ((II==JJ)-y{II})./rho.*g{1} + y{II}.*dg_du{1,JJ} );
            end
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

%Rusanov-ALE flux for convective term (from Badwaik et al. 2019, Mishra et al.):
%   f       = 0.5*(fL+fR) + 0.5*lambda*(uL-uR),
%   lambda  = max(|vL-wL|, |vR-wR|)
%In order to compute the Jacobian in an easier way, we write lambda as
%   lambda  = CL*SL*(vL-wL) + CR*SR*(vR-wR)
%Since fL=(vL-wL)*uL, fR=(vR-wR)*uR, then
%   f       = 0.5*(vL-wL)*uL + 0.5*(vR-wR)*uR +
%               0.5*(CL*SL*(vL-wL)+CR*SR*(vR-wR))* (uL-uR)
function [f, df_duL, df_duR] = f_Rusanov(model, uL, uR, ComputeJ)

    %Relative velocities:
    vwL     = uL{end-1}-uL{end};
    vwR     = uR{end-1}-uR{end};
    if abs(vwL)>=abs(vwR)
        CL  = 1.0;
        CR  = 0.0;
    else
        CL  = 0.0;
        CR  = 1.0;
    end
    SL      = sign(vwL);
    SR      = sign(vwR);
    lambda  = CL.*SL.*vwL+CR.*SR.*vwR;
    
    %Compute fluxes:
    f       = cell(model.nDiff, 1);
    df_duL  = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, uL{1});
    df_duR  = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, uR{1});
    for II=1:model.nDiff
        f{II}                   = 0.5.*vwL.*uL{II} + 0.5.*vwR.*uR{II} + ...
                                    0.5.*lambda.*(uL{II}-uR{II});
        if ComputeJ
            df_duL{II,II}       = + 0.5.*vwL + 0.5.*lambda;
            df_duL{II,end-1}    = + 0.5.*uL{II} + 0.5.*CL.*SL.*(uL{II}-uR{II});
            df_duL{II,end}      = - 0.5.*uL{II} - 0.5.*CL.*SL.*(uL{II}-uR{II});
            df_duR{II,II}       = + 0.5.*vwR - 0.5.*lambda;
            df_duR{II,end-1}    = + 0.5.*uR{II} + 0.5.*CR.*SR.*(uL{II}-uR{II});
            df_duR{II,end}      = - 0.5.*uR{II} - 0.5.*CR.*SR.*(uL{II}-uR{II});
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
    function Dmax=Dmaxfun(u)
        %its value is not necessary
        [rhoy, ~, ~, ~, y]  = aux_rho(model, u, du_dx);
        H                   = u{model.nSpecies+1};
        T                   = calc_T(H, rhoy, model);
        D_T                 = calc_D_T(T,y,model);
        D_rho               = calc_D_rho(T,y,model);
        Dmax                = max([D_T, D_rho]);
    end
    DmaxL                   = Dmaxfun(uL);
    DmaxR                   = Dmaxfun(uR);
    sigma                   = model.CW*max([DmaxL, DmaxR]);
    
    %Penalty terms:
    f       = cell(model.nDiff, 1);
    df_duL  = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, uL{1});
    df_duR  = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, uR{1});
    for II=1:model.nDiff
        f{II}               = sigma*(uL{II}-uR{II})./hp;
        if ComputeJ
            df_duL{II,II}   = + sigma./hp;
            df_duR{II,II}   = - sigma./hp;
        end
    end
    
end

function [f, df_duL, df_duL_dx, df_duR, df_duR_dx] = ...
    ftilde(model, ~, ~, uL, duL_dx, uR, duR_dx, hp, ComputeJ)
    
    %Evaluate convective, diffusive and penalty fluxes:
    [fc, dfc_duL, dfc_duR]      = f_Rusanov(model, uL, uR, ComputeJ);
    [fL, dfL_duL, dfL_duL_dx]   = f_diffusive(model, uL, duL_dx, ComputeJ);
    [fR, dfR_duR, dfR_duR_dx]   = f_diffusive(model, uR, duR_dx, ComputeJ);
    [fp, dfp_duL, dfp_duR]      = f_penalty(model, uL, uR, hp, ComputeJ);
    
    %Compute total flux:
    f           = cell(model.nDiff, 1);
    df_duL      = cell(model.nDiff, model.nVars);
    df_duL_dx   = cell(model.nDiff, model.nVars);
    df_duR      = cell(model.nDiff, model.nVars);
    df_duR_dx   = cell(model.nDiff, model.nVars);
    for II=1:model.nDiff
        f{II}   = fc{II} + 0.5*(fL{II} + fR{II}) + fp{II};
    end
    if ComputeJ
       for II=1:model.nDiff
           for JJ=1:model.nVars
               df_duL{II,JJ}    = dfc_duL{II,JJ} + 0.5*dfL_duL{II,JJ} + dfp_duL{II,JJ};
               df_duR{II,JJ}    = dfc_duR{II,JJ} + 0.5*dfR_duR{II,JJ} + dfp_duR{II,JJ};
               df_duL_dx{II,JJ} = 0.5*dfL_duL_dx{II,JJ};
               df_duR_dx{II,JJ} = 0.5*dfR_duR_dx{II,JJ};
           end
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
    u1              = model.u1(t);                      
    uL              = cell(model.nVars,1);
    for II=1:model.nDiff
        uL{II}      = u1{II};
    end
    uL{model.nDiff+1}   = u1{model.nDiff+1};
    uL{model.nVars}     = u{model.nVars};
    
    %Right state is the numerical solution:
    uR              = u;
    duR_dx          = du_dx; 
    
    %Convection based on left state:
    fc              = cell(model.nDiff,1);
    dfc_duL         = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, uL{1});
    dfc_duR         = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, uR{1});
    for II=1:model.nDiff
        fc{II}                  = (uL{end-1}-uL{end}).*uL{II};
        if ComputeJ
            dfc_duL{II,II}      = uL{end-1}-uL{end};
            dfc_duL{II,end-1}   = uL{II};
            dfc_duL{II,end}     = -uL{II};
        end
    end
    
    %Evaluate convective, diffusive and penalty fluxes:
%     [fc, dfc_duL, dfc_duR]      = f_Rusanov(model, uL, uR, ComputeJ);
    [fR, dfR_duR, dfR_duR_dx]   = f_diffusive(model, uR, duR_dx, ComputeJ);
    [fp, dfp_duL, dfp_duR]      = f_penalty(model, uL, uR, hp, ComputeJ);
    
%     disp('gas-flux1')
%     disp(fc)
%     disp(fR)
%     disp(fp)
    
    %Compute total flux:
    f           = cell(model.nDiff, 1);
    df_duL      = cell(model.nDiff, model.nVars);
    %df_duL_dx  = 0.0
    df_duR      = cell(model.nDiff, model.nVars);
    df_duR_dx   = cell(model.nDiff, model.nVars);
    for II=1:model.nDiff
        f{II}   = fc{II} + fR{II} + fp{II};
    end
    if ComputeJ
       for II=1:model.nDiff
           for JJ=1:model.nVars
               df_duL{II,JJ}    = dfc_duL{II,JJ} + dfp_duL{II,JJ};
               df_duR{II,JJ}    = dfc_duR{II,JJ} + dfR_duR{II,JJ} + dfp_duR{II,JJ};
               df_duR_dx{II,JJ} = dfR_duR_dx{II,JJ};
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
    ftildeN_Dirichlet(model, t, x, u, du_dx, hp, ComputeJ)
    
    %Left state is the numerical solution:
    uL              = u;
    duL_dx          = du_dx;    
    
    %Right state:
    %   -Differential variables are imposed from boundary conditions.
    %   -Velocity is extrapolated
    %   -Mesh velocity is extrapolated from solution
    %Diffusive flux is extrapolated:
    uN              = model.uN(t);                      
    uR              = cell(model.nVars,1);
    for II=1:model.nDiff
        uR{II}      = uN{II};
    end
    uR{model.nDiff+1}   = u{model.nDiff+1};
    uR{model.nVars}     = u{model.nVars};
    
    %Evaluate convective, diffusive and penalty fluxes:
    [fc, dfc_duL, dfc_duR]      = f_Rusanov(model, uL, uR, ComputeJ);
    [fL, dfL_duL, dfL_duL_dx]   = f_diffusive(model, uL, duL_dx, ComputeJ);
    [fp, dfp_duL, dfp_duR]      = f_penalty(model, uL, uR, hp, ComputeJ);
    
    %Compute total flux:
    f           = cell(model.nDiff, 1);
    df_duL      = cell(model.nDiff, model.nVars);
    df_duL_dx   = cell(model.nDiff, model.nVars);
    df_duR      = cell(model.nDiff, model.nVars);
    %df_duR_dx  = 0.0
    for II=1:model.nDiff
        f{II}   = fc{II} + fL{II} + fp{II};
    end
    if ComputeJ
       for II=1:model.nDiff
           for JJ=1:model.nVars
               df_duL{II,JJ}    = dfc_duL{II,JJ} + dfL_duL{II,JJ} +  dfp_duL{II,JJ};
               df_duR{II,JJ}    = dfc_duR{II,JJ} + dfp_duR{II,JJ};
               df_duL_dx{II,JJ} = dfL_duL_dx{II,JJ};
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

%df_dq means derivatives of flux w.r.t. the parameters that define the
%boundary condition (rhobar, Hbar in this case)
function [f, df_du, df_du_dx, df_dq] = ...
    ftildeN_NoFlux(model, t, x, u, du_dx, hp, ComputeJ)
    
    f       = Cells_Allocate(model.nDiff, 1, true, u{1});
    df_du   = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, u{1});
    df_du_dx= Cells_Allocate(model.nDiff, model.nVars, ComputeJ, u{1});
    df_dq   = Cells_Allocate(model.nDiff, 0, ComputeJ, u{1});
    
end

