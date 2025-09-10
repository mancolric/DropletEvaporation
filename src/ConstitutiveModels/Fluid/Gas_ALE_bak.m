function model=Gas_ALE_bak(fuel_names, comp_inerts, frac_masG)

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
    model.Tmin        = 300; 
    model.Tmax        = 1700;
    model.P           = 101325;
    model.N_polyfit   = 5;
    model.matrix      = Polyfit_properties_pureCompounds(model.gota, model.comp_inerts, model.Tmin, model.Tmax, model.P, model.N_polyfit);
    
    %Mandatory fields:
    model.nDiff       = 1+model.nSpecies;  %Number of differential variables
    model.nAlg        = 1;                 %Number of algebraic variables
    model.nVars       = model.nSpecies+3;  %Total nb of variables, including mesh velocity
    model.fQg         = @fQg;              %Function to compute flux, source terms and restriction
    model.ftilde      = @ftilde;           %Function to compute numerical flux at the internal faces
    model.ftilde1     = @ftilde1;          %Function to compute numerical flux at face 1, i.e., impose bondary condition
    model.ftildeN     = @ftildeN;          %Function to compute numerical flux at last face, i.e., impose bondary condition
    model.f_diffusive = @f_diffusive;
    
end
   
%--------------------------------------------------------------------------
%The functions below are only seen through structure model:

%Convection-diffusion flux and source terms:
function [  f, df_du, df_du_dx, ...
            Q, dQ_du, dQ_du_dx, ...
            g, dg_du, dg_du_dx ] = ...
    fQg(model, t, x, u, du_dx, ComputeJ)

    %Extract variables:
    [rhoy, ~, rho, ~, y] = aux_rho(model, u, du_dx);
    H           = u{model.nDiff};
    v           = u{end-1};
    w           = u{end};
    
    %Diffusive flux: 
    [f, df_du, df_du_dx]    = f_diffusive(model, u, du_dx, ComputeJ);
    %Add convective terms:
    for II=1:model.nDiff
        f{II}               = f{II} + (v-w).*u{II};
        if ComputeJ
            df_du{II,II}    = df_du{II,II} + v-w;
            df_du{II,end-1} = df_du{II,end-1} + u{II};
            df_du{II,end}   = df_du{II,end} - u{II};
        end
    end
    
    T       = calc_T(H, rhoy, model);
    rho_bar = calc_rho(model, y, T);

    %Source and restriction:
    Q                       = model.Q(t,x);
    g                       = { rho - rho_bar };
    [dQ_du, dQ_du_dx]       = Cells_Allocate(model.nDiff, model.nVars, ComputeJ, rhoy{1});
    [dg_du, dg_du_dx]       = Cells_Allocate(model.nAlg, model.nVars, ComputeJ, rhoy{1});
    if ComputeJ
        
        %Derivatives w.r.t. rho and H:
        for II=1:model.nSpecies
            dg_du{1,II}         = 1.0 + 0.0*u{1}; 
        end
        
    end
end

function [rhoy, drhoy_dx, rho, sum_drhoy_dx, y] = aux_rho(model, u, du_dx)

    %Extract variables:
    rhoy        = cell(model.nSpecies, 1);
    drhoy_dx    = cell(model.nSpecies, 1);
    for II=1:model.nSpecies
        rhoy{II}     = u{II};
        drhoy_dx{II} = du_dx{II};
    end
   
    rho          = zeros(size(rhoy{1})); % Mezcla
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

function [f, df_du, df_du_dx] = ...
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

end

function [f, df_duL, df_duL_dx, df_duR, df_duR_dx] = ...
    ftilde(model, ~, ~, uL, duL_dx, uR, duR_dx, hp, ComputeJ)
    
    %Allocate derivatives:
    df_duL      = cell(model.nDiff, model.nVars);
    df_duL_dx   = cell(model.nDiff, model.nVars);
    df_duR      = cell(model.nDiff, model.nVars);
    df_duR_dx   = cell(model.nDiff, model.nVars);
    
    %Evaluate diffusive flux f_diff(0.5*(uL+uR), 0.5*(duL_dx+duR_dx)):
    u                       = cell(model.nVars, 1);
    du_dx                   = cell(model.nVars, 1);
    for II=1:model.nVars
        u{II}               = 0.5*(uL{II} + uR{II});
        du_dx{II}           = 0.5*(duL_dx{II} + duR_dx{II});
    end
    [f, df_du, df_du_dx]    = f_diffusive(model, u, du_dx, ComputeJ);
    if ComputeJ
        
        %Apply chain rule:
        for II=1:model.nDiff
            for JJ=1:model.nVars
                df_duL{II,JJ}       = df_du{II,JJ}*0.5;
                df_duR{II,JJ}       = df_du{II,JJ}*0.5;
                df_duL_dx{II,JJ}    = df_du_dx{II,JJ}*0.5;
                df_duR_dx{II,JJ}    = df_du_dx{II,JJ}*0.5;
            end
        end
        
    end
    
    %Diffusion due to penalty:
    [rhoy, ~, ~, ~, y]   = aux_rho(model, u, du_dx);
    H           = u{model.nSpecies+1};
    T           = calc_T(H, rhoy, model);
    D_T         = calc_D_T(T,y,model);
    D_rho       = calc_D_rho(T,y,model);
    sigma       = model.CW*max([D_T; D_rho]);
    
    %Rusanov-ALE flux for convective term (from Badwaik et al. 2019, Mishra et al.):
    %   f       = 0.5*(fL+fR) + 0.5*lambda*(uL-uR),
    %   lambda  = max(|vL-wL|, |vR-wR|)
    %In order to compute the Jacobian in an easier way, we write lambda as
    %   lambda  = CL*SL*(vL-wL) + CR*SR*(vR-wR)
    %Since fL=(vL-wL)*uL, fR=(vR-wR)*uR, then
    %   f       = 0.5*(vL-wL)*uL + 0.5*(vR-wR)*uR +
    %               0.5*(CL*SL*(vL-wL)+CR*SR*(vR-wR))* (uL-uR)
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
    for II=1:model.nDiff
        f{II}                   = f{II} + 0.5.*vwL.*uL{II} + 0.5.*vwR.*uR{II} + ...
                                    0.5.*lambda.*(uL{II}-uR{II});
        if ComputeJ
            df_duL{II,II}       = df_duL{II,II} + 0.5.*vwL + 0.5.*lambda;
            df_duL{II,end-1}    = df_duL{II,end-1} + 0.5.*uL{II} + 0.5.*CL.*SL.*(uL{II}-uR{II});
            df_duL{II,end}      = df_duL{II,end} - 0.5.*uL{II} - 0.5.*CL.*SL.*(uL{II}-uR{II});
            df_duR{II,II}       = df_duR{II,II} + 0.5.*vwR - 0.5.*lambda;
            df_duR{II,end-1}    = df_duR{II,end-1} + 0.5.*uR{II} + 0.5.*CR.*SR.*(uL{II}-uR{II});
            df_duR{II,end}      = df_duR{II,end} - 0.5.*uR{II} - 0.5.*CR.*SR.*(uL{II}-uR{II});
        end
    end
    
    %Add penalty terms:
    for II=1:model.nDiff
        f{II}               = f{II} + sigma*(uL{II}-uR{II})./hp;
        if ComputeJ
            df_duL{II,II}   = df_duL{II,II} + sigma./hp;
            df_duR{II,II}   = df_duR{II,II} - sigma./hp;
        end
    end
        
end

%df_dq means derivatives of flux w.r.t. the parameters that define the
%boundary condition (rhobar, Hbar, vbar in this case)
function [f, df_du, df_du_dx, df_dq] = ...
    ftilde1(model, t, x, u, du_dx, hp, ComputeJ)
    
    %Left state:
    %   -Differential variables are imposed from boundary conditions.
    %   -Velocity is imposed
    %   -Mesh velocity is extrapolated from solution
    %   -Derivatives are extrapolated
    u1              = model.u1(t);                      
    uL              = cell(model.nVars,1);
    for II=1:model.nDiff
        uL{II}      = u1{II};
    end
    uL{model.nDiff+1}   = u1{model.nDiff+1};
    uL{model.nVars}     = u{model.nVars};
    duL_dx              = du_dx; 
    
    %Right state is the numerical solution:
    uR              = u;
    duR_dx          = du_dx;    
    
    %Evaluate numerical flux:
    [f, df_duL, df_duL_dx, df_duR, df_duR_dx]  = ftilde(model, t, x, uL, duL_dx, uR, duR_dx, hp, ComputeJ);
    df_du           = cell(model.nDiff, model.nVars);
    df_du_dx        = cell(model.nDiff, model.nVars);
    df_dq           = cell(model.nDiff, model.nDiff+model.nAlg);  %Derivatives w.r.t. rhobarY_i, Hbar, vbar
    if ComputeJ
        for II=1:model.nDiff
            for JJ=1:model.nVars
                df_du{II,JJ}    = df_duL{II,JJ}*(JJ==model.nVars) + df_duR{II,JJ}*1.0;
                df_du_dx{II,JJ} = df_duL_dx{II,JJ}*1.0 + df_duR_dx{II,JJ}*1.0;
            end
            for JJ=1:model.nDiff+model.nAlg
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
    %   -Derivatives are extrapolated
    uN              = model.uN(t);                      
    uR              = cell(model.nVars,1);
    for II=1:model.nDiff
        uR{II}      = uN{II};
    end
    uR{model.nDiff+1}   = u{model.nDiff+1};
    uR{model.nVars}     = u{model.nVars};
    duR_dx              = du_dx; 
    
    %Evaluate numerical flux:
    [f, df_duL, df_duL_dx, df_duR, df_duR_dx]  = ftilde(model, t, x, uL, duL_dx, uR, duR_dx, hp, ComputeJ);
    df_du           = cell(model.nDiff, model.nVars);
    df_du_dx        = cell(model.nDiff, model.nVars);
    df_dq           = cell(model.nDiff, model.nDiff);  %Derivatives w.r.t. rhobarY_i, Hbar
    if ComputeJ
        for II=1:model.nDiff
            for JJ=1:model.nVars
                df_du{II,JJ}    = df_duL{II,JJ}*1.0 + df_duR{II,JJ}*(JJ==model.nVars-1) + df_duR{II,JJ}*(JJ==model.nVars);
                df_du_dx{II,JJ} = df_duL_dx{II,JJ}*1.0 + df_duR_dx{II,JJ}*1.0;
            end
            for JJ=1:model.nDiff
                df_dq{II,JJ}     = df_duR{II,JJ};
            end
        end
    end
    
end
