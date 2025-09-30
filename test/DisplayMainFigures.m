format long
clear
clc
tic
close all

%LOAD RESULTS:
aux1 = load('Octano0.7+Eicosano.mat');

%CALCULATE EVAPORATED MASS FLOW:
m_l_1 = 4*pi().*((cell2mat(aux1.d)./2).^2).*cell2mat(aux1.m_l);     

% CELL2MAT t:
t_n = cell2mat(aux1.t_n);

%PLOT EVAPORATED MASS FLOW:
subplot(1,3,1)
hold off
plot(t_n, m_l_1, 'r--', 'LineWidth', 1.5)
hold on
set(gca, 'FontSize', 15)
ylabel('$\dot{m}$ [kg/s]', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold')
xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold')
% legend('FontSize', 12, 'Location', 'Best', 'Interpreter', 'latex')

%PLOT d^2:
subplot(1,3,2)
hold off
plot(t_n, cell2mat(aux1.d).^2, 'r--', 'LineWidth', 2)
hold on
set(gca, 'FontSize', 15)
ylabel('$d^2 \, [\mathrm{m}^2]$', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold')
xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold')
% legend('FontSize', 12, 'Location', 'Best')

%PLOT SURFACE TEMPERATURE:
subplot(1,3,3)
hold off
plot(t_n, cell2mat(aux1.T_ql), 'r--', 'LineWidth', 1.5)
hold on
set(gca, 'FontSize', 15)
ylabel('$T_{\mathrm{s}}$ [K]', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold')
xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold')
% legend('FontSize', 12, 'Location', 'Best')