clear
clc
close all
load('./results/Saved_varsEnd.mat')

Nt              = length([save_vars{6,:}]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%              m_dot                 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mdot_i          = [save_vars{6,1:Nt}];
mdot_total      = sum(cell2mat(mdot_i),1);

figure(1)
plot([save_vars{1, 1:Nt}], mdot_total, '-b')
hold on

xlabel("$$t[s]$$", Interpreter="latex", FontSize=16)
ylabel("$$\dot{m} [Kg/m^2s]$$", Interpreter="latex", FontSize=16)
legend("Present model", 'Interpreter','latex', 'Location','northwest', FontSize=16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%              T_gas                 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time_instant    = 3;
disp("Time instant " + Time_instant + " correspond to t = " + save_vars{1,Time_instant})


T_gas           = [save_vars{12,Time_instant}];
x_gas           = [save_vars{17,Time_instant}];

figure(2)
plot(x_gas, T_gas, '-b')
hold on

xlabel("$$r[m]$$", Interpreter="latex", FontSize=16)
ylabel("$$Temperature [K]$$", Interpreter="latex", FontSize=16)
legend("Present model", 'Interpreter','latex', 'Location','northwest', FontSize=16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%              Y_liq                 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time_instant    = 16;
disp("Time instant " + Time_instant + " correspond to t = " + save_vars{1,Time_instant})

Y_liq           = save_vars{4,Time_instant};
nSpecies        = length(Y_liq);
names           = cell(nSpecies,1);
x_liq           = [save_vars{16,Time_instant}];

figure(3)
hold on
for ii=1:nSpecies
    plot(x_liq, Y_liq{ii})
    name        = "Specie " + ii;
    names{ii}   = name;
end

xlabel("$$r[m]$$", Interpreter="latex", FontSize=16)
ylabel("$$Temperature [K]$$", Interpreter="latex", FontSize=16)
legend(names, 'Interpreter','latex', 'Location','northwest', 'FontSize', 16);