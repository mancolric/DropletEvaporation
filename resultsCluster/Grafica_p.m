format long
clear
clc
tic
close all

isFilled = @(c) ~(isempty(c) || (isnumeric(c) && isscalar(c) && isnan(c)));
trim_all_empty_cols = @(C) C(:, any(cellfun(isFilled, C), 1));

aux1 = load('Saved_varsEnd(h).mat');
aux2 = load('Saved_varsEnd(h_2).mat');
aux3 = load('Saved_varsEnd(h_4).mat');
% aux13 = load('Saved_varsEnd(h_8).mat');
aux1.save_vars = trim_all_empty_cols(aux1.save_vars);
aux2.save_vars = trim_all_empty_cols(aux2.save_vars);
aux3.save_vars = trim_all_empty_cols(aux3.save_vars);

aux4 = load('Saved_varsEnd2(h).mat');
aux5 = load('Saved_varsEnd2(h_2).mat');
aux6 = load('Saved_varsEnd2(h_4).mat');
% aux14 = load('Saved_varsEnd2(h_8).mat');
aux4.save_vars = trim_all_empty_cols(aux4.save_vars);
aux5.save_vars = trim_all_empty_cols(aux5.save_vars);
aux6.save_vars = trim_all_empty_cols(aux6.save_vars);

aux7 = load('Saved_varsEnd3(h).mat');
aux8 = load('Saved_varsEnd3(h_2).mat');
aux9 = load('Saved_varsEnd3(h_4).mat');
aux7.save_vars = trim_all_empty_cols(aux7.save_vars);
aux8.save_vars = trim_all_empty_cols(aux8.save_vars);
aux9.save_vars = trim_all_empty_cols(aux9.save_vars);

% aux10 = load('Saved_varsEnd4(h).mat');
% aux11 = load('Saved_varsEnd4(h_2).mat');
% aux12 = load('Saved_varsEnd4(h_4).mat');

row = 2;

% === GRUPO 1 ===
T_h   = aux1.save_vars{row, end};
T_h_2 = aux2.save_vars{row, end};
T_h_4 = aux3.save_vars{row, end};
% T_h_8 = aux13.save_vars{row, 1912};

alpha1 = log2(abs(T_h - T_h_2) ./ abs(T_h_2 - T_h_4));
T_exact1 = (T_h - (2.^alpha1).*T_h_2) ./ (1 - (2.^alpha1));

err1_h   = abs(T_h   - T_exact1)/T_exact1;
err1_h2  = abs(T_h_2 - T_exact1)/T_exact1;
err1_h4  = abs(T_h_4 - T_exact1)/T_exact1;
% err1_h8  = abs(T_h_8 - T_exact1)/T_exact1;
% err1 = [err1_h, err1_h2, err1_h4, err1_h8];
err1 = [err1_h, err1_h2, err1_h4];

tCPU1_h   = aux1.simulationTime;
tCPU1_h2  = aux2.simulationTime;
tCPU1_h4  = aux3.simulationTime;
% tCPU1_h8  = aux13.simulationTime;

% tCPU1 = [tCPU1_h,tCPU1_h2,tCPU1_h4,tCPU1_h8];
tCPU1 = [tCPU1_h,tCPU1_h2,tCPU1_h4];

% === GRUPO 2 ===
T_h2   = aux4.save_vars{row, end};
T_h2_2 = aux5.save_vars{row, end};
T_h2_4 = aux6.save_vars{row, end};
% T_h2_8 = aux14.save_vars{row, 4849};

alpha2 = log2(abs(T_h2 - T_h2_2) / abs(T_h2_2 - T_h2_4));
T_exact2 = (T_h2- (2.^alpha2).*T_h2_2) ./ (1 - (2.^alpha2));

err2_h   = abs(T_h2   - T_exact2)/T_exact2;
err2_h2  = abs(T_h2_2 - T_exact2)/T_exact2;
err2_h4  = abs(T_h2_4 - T_exact2)/T_exact2;
% err2_h8  = abs(T_h_8 - T_exact2)/T_exact2;
% err2 = [err2_h, err2_h2, err2_h4, err2_h8];
err2 = [err2_h, err2_h2, err2_h4];

tCPU2_h   = aux4.simulationTime;
tCPU2_h2  = aux5.simulationTime;
tCPU2_h4  = aux6.simulationTime;
% tCPU2_h8  = aux14.simulationTime;

% tCPU2 = [tCPU2_h,tCPU2_h2,tCPU2_h4,tCPU2_h8];
tCPU2 = [tCPU2_h,tCPU2_h2,tCPU2_h4];

% === GRUPO 3 ===
T_h3   = aux7.save_vars{row, end};
T_h3_2 = aux8.save_vars{row, end};
T_h3_4 = aux9.save_vars{row, end};

alpha3 = log2(abs(T_h3 - T_h3_2) ./ abs(T_h3_2 - T_h3_4));
T_exact3 = (T_h3 - (2.^alpha3).*T_h3_2) ./ (1 - (2.^alpha3));

err3_h   = abs(T_h3   - T_exact3)/T_exact3;
err3_h2  = abs(T_h3_2 - T_exact3)/T_exact3;
err3_h4  = abs(T_h3_4 - T_exact3)/T_exact3;
err3 = [err3_h, err3_h2, err3_h4];

tCPU3_h   = aux7.simulationTime;
tCPU3_h2  = aux8.simulationTime;
tCPU3_h4  = aux9.simulationTime;

tCPU3 = [tCPU3_h,tCPU3_h2,tCPU3_h4];

% % === GRUPO 4 ===
% T_h   = aux10.save_vars{row, };
% T_h_2 = aux11.save_vars{row, };
% T_h_4 = aux12.save_vars{row, };
% 
% alpha1 = log2(abs(T_h - T_h_2) ./ abs(T_h_2 - T_h_4));
% T_exact4 = (T_h - (2.^alpha1).*T_h_2) ./ (1 - (2.^alpha1));
% 
% err4_h   = abs(T_h   - T_exact1)/T_exact1;
% err4_h2  = abs(T_h_2 - T_exact1)/T_exact1;
% err4_h4  = abs(T_h_4 - T_exact1)/T_exact1;
% err4 = [err4_h, err4_h2, err4_h4];
% 
% tCPU4 = [,,];

% === GRÁFICA ===
figure;
hold on; grid on; box on;

plot(tCPU1, err1, '-s', 'MarkerFaceColor', 'b', 'LineWidth', 2, 'DisplayName', 'p = 2', 'MarkerEdgeColor', 'b', 'Color', 'b');
plot(tCPU2, err2, '-s', 'MarkerFaceColor', 'r', 'LineWidth', 2, 'DisplayName', 'p = 3', 'MarkerEdgeColor', 'r', 'Color', 'r');
plot(tCPU3, err3, '-s', 'MarkerFaceColor', 'g', 'LineWidth', 2, 'DisplayName', 'p = 4', 'MarkerEdgeColor', 'g', 'Color', 'g');
% plot(tCPU4, err4, '-s', 'MarkerFaceColor', 'y', 'LineWidth', 2, 'DisplayName', 'p = 5', 'MarkerEdgeColor', 'y', 'Color', 'y');

xlabel('$t_{CPU}\,[\mathrm{s}]$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$e_{ST}$', 'Interpreter', 'latex', 'FontSize', 16);
% title('Error frente al tiempo de simulación', 'FontSize', 16);

legend('Location', 'best', 'FontSize', 14);

set(gca, 'XScale', 'log', 'YScale', 'log');  % escalas logarítmicas en ambos ejes
set(gca, 'FontSize', 14);  % agranda los números de los ejes
grid on;                % activa la cuadrícula principal
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');  % desactiva las cuadrículas menores
box on;                 % mantiene el borde del gráfico

% === RESULTADOS ===
Resultados = table( ...
    [T_h;   T_h2;   T_h3], ...
    [T_h_2; T_h2_2; T_h3_2], ...
    [T_h_4; T_h2_4; T_h3_4], ...
    [T_exact1; T_exact2; T_exact3], ...
    [alpha1; alpha2; alpha3], ...
    'VariableNames', {'T_h', 'T_h_2', 'T_h_4', 'T_exact', 'alpha'}, ...
    'RowNames', {'Caso 1', 'Caso 2', 'Caso 3'});

disp('-------------------------------------------');
disp('        RESULTADOS DE TEMPERATURAS Y α     ');
disp('-------------------------------------------');
disp(Resultados);