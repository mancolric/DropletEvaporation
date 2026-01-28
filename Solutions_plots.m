clear
clc
close all
load('./results/Saved_varsEnd.mat')

mdot_Alvaro = [save_vars{6, 1:24}];
mdot_gas_D  = [save_vars{17,1:24}];
mdot_gas_D  = sum(cell2mat(mdot_gas_D),1);
mdot_liq_D  = [save_vars{22,1:24}];
mdot_liq_D  = sum(cell2mat(mdot_liq_D),1);

figure(1)
plot([save_vars{1, 1:24}], mdot_Alvaro, '-b')
hold on
plot([save_vars{1, 1:24}], mdot_gas_D, '-r')
plot([save_vars{1, 1:24}], mdot_liq_D, '.-g')

xlabel("$$t[s]$$", Interpreter="latex", FontSize=16)
ylabel("$$\dot{m} [Kg/m^2s]$$", Interpreter="latex", FontSize=16)
legend("$$\dot{m}$$ origianl", "$$\dot{m}$$ nuevo desde gas", "$$\dot{m}$$ nuevo desde liquido", 'Interpreter','latex', 'Location','northwest', FontSize=16)