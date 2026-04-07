clear
clc
close all
fuel_names          = {'Heptano'};
inert_comps         = {'N2', 'O2', 'CO2', 'H2O'};
mass_fracG          = {0.79, 0.21, 0.0, 0.0};
addpath(genpath(pwd));


model               = Gas_ALE(fuel_names, inert_comps, mass_fracG);

T = 1000;
rhoy = {0,0.8,0,0,0.2}';

Q = Qfun(model, T, rhoy);
for ii=1:length(Q)
    disp(Q{ii})
end

T = ones(4)*1000;
rhoy = {ones(4)*0,ones(4)*0.8,ones(4)*0,ones(4)*0,ones(4)*0.2}';

Q2 = Qfun(model, T, rhoy);
for ii=1:length(Q2)
    disp(Q2{ii}(1,1)-Q{ii})
end