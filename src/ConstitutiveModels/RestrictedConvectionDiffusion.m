%Create linear convection-diffusion model. Some parameters are input and
%some others are default.

function model = RestrictedConvectionDiffusion(epsilon)

    %User-defined parameters:
    model.epsilon   = epsilon;              %Viscosity
    model.CW        = 50.0;                 %Weight coefficient for penalty terms
    model.Q         = @(t,x) {0.0*x};       %Source term that depends on (t,x)
    model.u1        = @(t,x) {2.0; 1.0};    %Boundary condition at the left
    model.uN        = @(t,x) {1.0; NaN};    %Boundary condition at the right
    model.u1target  = @(t,x) 1.0./x;        %Restriction for u1
    
    %Mandatory fields:
    model.nDiff     = 1;            %Number of differential variables
    model.nAlg      = 1;            %Number of algebraic variables
    model.fQg       = @fQg;         %Function to compute flux, source terms and restriction
    model.ftilde    = @ftilde;      %Function to compute numerical flux at the internal faces
    model.ftilde1   = @ftilde1;     %Function to compute numerical flux at face 1, i.e., impose bondary condition
    model.ftildeN   = @ftildeN;     %Function to compute numerical flux at last face, i.e., impose bondary condition
    
end

%--------------------------------------------------------------------------
%The functions below are only seen through structure model:

%Convection-diffusion flux and source terms:
function [  f, df_du, df_du_dx, ...
            Q, dQ_du, dQ_du_dx, ...
            g, dg_du, dg_du_dx ] = ...
    fQg(model, t, x, u, du_dx, ComputeJ)

    %Flux, source and restriction:
    f           = { u{2}.*u{1} - model.epsilon*du_dx{1} };
    Q           = model.Q(t,x);
    g           = { u{1} - model.u1target(t,x) };
    
    %Derivatives:
    df_du       = cell(model.nDiff,model.nDiff+model.nAlg);
    df_du_dx    = cell(model.nDiff,model.nDiff+model.nAlg);
    dQ_du       = cell(model.nDiff,model.nDiff+model.nAlg);
    dQ_du_dx    = cell(model.nDiff,model.nDiff+model.nAlg);
    dg_du       = cell(model.nAlg,model.nDiff+model.nAlg);
    dg_du_dx    = cell(model.nAlg,model.nDiff+model.nAlg);
    if ComputeJ
        %Derivatives w.r.t. u1:
        df_du{1,1}      = u{2};
        df_du_dx{1,1}   = 0.0*u{1} - model.epsilon;
        dQ_du{1,1}      = 0.0*u{1};
        dQ_du_dx{1,1}   = 0.0*u{1};
        dg_du{1,1}      = 0.0*u{1} + 1.0;
        dg_du_dx{1,1}   = 0.0*u{1};
        
        %Derivatives w.r.t. u2:
        df_du{1,2}      = u{1};
        df_du_dx{1,2}   = 0.0*u{1};
        dQ_du{1,2}      = 0.0*u{1};
        dQ_du_dx{1,2}   = 0.0*u{1};
        dg_du{1,2}      = 0.0*u{1};
        dg_du_dx{1,2}   = 0.0*u{1};
    end
end

function [f, df_duL, df_duL_dx, df_duR, df_duR_dx] = ...
    ftilde(model, t, x, uL, duL_dx, uR, duR_dx, h, ComputeJ)
    
    %ftilde = f(u, du/dx) + CW*epsilon*(uL-uR), with
    %   u=uL, du/dx = 0.5*(duL/dx + duR/dx)
    u           = uL;
    du_dx       = cell(model.nDiff+model.nAlg);
    for II=1:model.nDiff+model.nAlg
        du_dx{II}   = 0.5*(duL_dx{II} + duR_dx{II});
    end
    
    %Evaluate flux and add penalty terms:
    [f, df_du, df_du_dx]    = fQg(model, t, x, u, du_dx, ComputeJ);
    f{1}                    = f{1} + model.CW*model.epsilon*(uL{1}-uR{1})./h;
                
    %Derivatives:
    df_duL      = cell(model.nDiff, model.nDiff+model.nAlg);
    df_duL_dx   = cell(model.nDiff, model.nDiff+model.nAlg);
    df_duR      = cell(model.nDiff, model.nDiff+model.nAlg);
    df_duR_dx   = cell(model.nDiff, model.nDiff+model.nAlg);
    if ComputeJ
        
        %Apply chain rule:
        for II=1:model.nDiff
            for JJ=1:model.nDiff+model.nAlg
                df_duL{II,JJ}       = df_du{II,JJ}*1.0 + model.CW*model.epsilon*(II==JJ)./h;
                df_duR{II,JJ}       = df_du{II,JJ}*0.0 - model.CW*model.epsilon*(II==JJ)./h;
                df_duL_dx{II,JJ}    = df_du_dx{II,JJ}*0.5;
                df_duR_dx{II,JJ}    = df_du_dx{II,JJ}*0.5;
            end
        end
        
    end
        
end

function [f, df_du, df_du_dx] = ...
    ftilde1(model, t, x, u, du_dx, h, ComputeJ)
    
    %Left and right states:
    uL          = model.u1(t,x);        %Impose Dirichlet condition
    uR          = u;
    duL_dx      = du_dx;                %Extrapolate du_dx
    duR_dx      = du_dx;    
    
    %Evaluate numerical flux:
    [f, df_duL, df_duL_dx, df_duR, df_duR_dx]  = ftilde(model, t, x, uL, duL_dx, uR, duR_dx, h, ComputeJ);
    df_du       = cell(model.nDiff, model.nDiff+model.nAlg);
    df_du_dx    = cell(model.nDiff, model.nDiff+model.nAlg);
    if ComputeJ
        for II=1:model.nDiff
            for JJ=1:model.nDiff+model.nAlg
                df_du{II,JJ}    = df_duR{II,JJ};
                df_du_dx{II,JJ} = df_duL_dx{II,JJ} + df_duR_dx{II,JJ};
            end
        end
    end
    
end

function [f, df_du, df_du_dx] = ...
    ftildeN(model, t, x, u, du_dx, h, ComputeJ)
    
    %Left and right states:
    uL          = u;
    uR          = model.uN(t,x);        %Impose Dirichlet condition
    duL_dx      = du_dx;                
    duR_dx      = du_dx;                %Extrapolate du_dx
    
    %Evaluate numerical flux:
    [f, df_duL, df_duL_dx, df_duR, df_duR_dx]  = ftilde(model, t, x, uL, duL_dx, uR, duR_dx, h, ComputeJ);
    df_du       = cell(model.nDiff, model.nDiff+model.nAlg);
    df_du_dx    = cell(model.nDiff, model.nDiff+model.nAlg);
    if ComputeJ
        for II=1:model.nDiff
            for JJ=1:model.nDiff+model.nAlg
                df_du{II,JJ}    = df_duL{II,JJ};
                df_du_dx{II,JJ} = df_duL_dx{II,JJ} + df_duR_dx{II,JJ};
            end
        end
    end
    
end
