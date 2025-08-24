%Create linear convection-diffusion model. Some parameters are input and
%some others are default.

function model = LinearConvectionDiffusion(a, epsilon)

    %User-defined parameters:
    model.a         = a;                %Velocity
    model.epsilon   = epsilon;          %Viscosity
    model.CW        = 50.0;             %Weight coefficient for penalty terms
    model.Q         = @(t,x) {0.0*x};   %Source term that depends on (t,x)
    model.u1        = @(t,x) {2.0};     %Boundary condition at the left
    model.uN        = @(t,x) {1.0};     %Boundary condition at the right
    
    %Mandatory fields:
    model.nVars     = 1;
    model.fQ        = @fQ;          %Function to compute flux and source terms
    model.ftilde    = @ftilde;      %Function to compute numerical flux at the internal faces
    model.ftilde1   = @ftilde1;     %Function to compute numerical flux at face 1, i.e., impose bondary condition
    model.ftildeN   = @ftildeN;     %Function to compute numerical flux at last face, i.e., impose bondary condition
    
end

%--------------------------------------------------------------------------
%The functions below are only seen through structure model:

%Convection-diffusion flux and source terms:
function [f, df_du, df_du_dx, Q, dQ_du, dQ_du_dx] = ...
    fQ(model, t, x, u, du_dx, ComputeJ)

    f           = { model.a*u{1} - model.epsilon*du_dx{1} };
    df_du       = cell(1,1);
    df_du_dx    = cell(1,1);
    Q           = model.Q(t,x);
    dQ_du       = cell(1,1);
    dQ_du_dx    = cell(1,1);
    if ComputeJ
        df_du{1,1}      = repmat(model.a, size(u{1}));
        df_du_dx{1,1}   = repmat(-model.epsilon, size(u{1}));
        dQ_du{1,1}      = zeros(size(u{1}));
        dQ_du_dx{1,1}   = zeros(size(u{1}));
    end
end

function [f, df_duL, df_duL_dx, df_duR, df_duR_dx] = ...
    ftilde(model, t, x, uL, duL_dx, uR, duR_dx, h, ComputeJ)
    
    %Upwind for convective term + average of diffusive terms:
    f           = { model.a*uL{1} - ...
                    model.epsilon*0.5*(duL_dx{1}+duR_dx{1}) + ...
                    model.CW*model.epsilon*(uL{1}-uR{1})./h };
    df_duL      = cell(1,1);
    df_duL_dx   = cell(1,1);
    df_duR      = cell(1,1);
    df_duR_dx   = cell(1,1);
    if ComputeJ
        df_duL{1,1}     = model.a + model.CW*model.epsilon./h;
        df_duL_dx{1,1}  = repmat(-0.5*model.epsilon, size(uL{1}));
        df_duR{1,1}     = -model.CW*model.epsilon./h;
        df_duR_dx{1,1}  = repmat(-0.5*model.epsilon, size(uR{1}));
    end
    
end

function [f, df_du, df_du_dx] = ...
    ftilde1(model, t, x, u, du_dx, h, ComputeJ)
    
    %Upwind for convective term + average of diffusive terms:
    uDir        = model.u1(t,x);
    f           = { model.a*uDir{1} - ...
                    model.epsilon*du_dx{1} + ...
                    model.CW*model.epsilon*(uDir{1}-u{1})./h };
    df_du       = cell(1,1);
    df_du_dx    = cell(1,1);
    if ComputeJ
        df_du{1,1}      = -model.CW*model.epsilon./h;
        df_du_dx{1,1}   = repmat(-model.epsilon, size(u{1}));
    end
    
end

function [f, df_du, df_du_dx] = ...
    ftildeN(model, t, x, u, du_dx, h, ComputeJ)
    
    %Upwind for convective term + average of diffusive terms:
    uDir        = model.uN(t,x);
    f           = { model.a*u{1} - ...
                    model.epsilon*du_dx{1} + ...
                    model.CW*model.epsilon*(u{1}-uDir{1})./h };
    df_du       = cell(1,1);
    df_du_dx    = cell(1,1);
    if ComputeJ
        df_du{1,1}      = model.a + model.CW*model.epsilon./h;
        df_du_dx{1,1}   = repmat(-model.epsilon, size(u{1}));
    end
    
end
