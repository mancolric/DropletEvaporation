format long
clear
clc
tic
close all

tq = 0.04145;

aux1 = load('Saved_varsEnd(h).mat');
aux2 = load('Saved_varsEnd(h_2).mat');
aux3 = load('Saved_varsEnd(h_4).mat');
% aux13 = load('Saved_varsEnd(h_8).mat');

aux4 = load('Saved_varsEnd2(h).mat');
aux5 = load('Saved_varsEnd2(h_2).mat');
aux6 = load('Saved_varsEnd2(h_4).mat');
% aux14 = load('Saved_varsEnd2(h_8).mat');

aux7 = load('Saved_varsEnd3(h).mat');
aux8 = load('Saved_varsEnd3(h_2).mat');
aux9 = load('Saved_varsEnd3(h_4).mat');

% aux10 = load('Saved_varsEnd4(h).mat');
% aux11 = load('Saved_varsEnd4(h_2).mat');
% aux12 = load('Saved_varsEnd4(h_4).mat');

% aux3.save_vars{2, 1547} = 5.05908843356480*1e2;
% save('Saved_varsEnd(h_4).mat', '-struct', 'aux3');

row = 2;

% === GRUPO 1 ===
T_h   = interp_from_aux(aux1, tq);
T_h_2 = interp_from_aux(aux2, tq);
T_h_4 = interp_from_aux(aux3, tq);

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
T_h   = interp_from_aux(aux4, tq);
T_h_2 = interp_from_aux(aux5, tq);
T_h_4 = interp_from_aux(aux6, tq);

alpha2 = log2(abs(T_h - T_h_2) ./ abs(T_h_2 - T_h_4));
T_exact2 = (T_h - (2.^alpha2).*T_h_2) ./ (1 - (2.^alpha2));

err2_h   = abs(T_h   - T_exact2)/T_exact2;
err2_h2  = abs(T_h_2 - T_exact2)/T_exact2;
err2_h4  = abs(T_h_4 - T_exact2)/T_exact2;
% err2_h8  = abs(T_h_8 - T_exact2)/T_exact2;
% err2 = [err2_h, err2_h2, err2_h4, err2_h8];
err2 = [err2_h, err2_h2, err2_h4];

tCPU2_h   = aux4.simulationTime;
tCPU2_h2  = aux5.simulationTime;
tCPU2_h4  = aux6.simulationTime;
% tCPU2_h8  = aux14.simulationTime;

% tCPU2 = [tCPU2_h,tCPU2_h2,tCPU2_h4,tCPU2_h8];
tCPU2 = [tCPU1_h,tCPU1_h2,tCPU1_h4];

% === GRUPO 3 ===
T_h   = interp_from_aux(aux7, tq);
T_h_2 = interp_from_aux(aux8, tq);
T_h_4 = interp_from_aux(aux9, tq);

alpha3 = log2((abs(T_h - T_h_2)) ./ abs(T_h_2 - T_h_4));
T_exact3 = (T_h - (2.^alpha3).*T_h_2) ./ (1 - (2.^alpha3));

err3_h   = abs(T_h   - T_exact3)/T_exact3;
err3_h2  = abs(T_h_2 - T_exact3)/T_exact3;
err3_h4  = abs(T_h_4 - T_exact3)/T_exact3;
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

% === RESULTADOS ===
disp('alpha1 ='); disp(alpha1);
disp('T_exact1 ='); disp(T_exact1);

disp('alpha2 ='); disp(alpha2);
disp('T_exact2 ='); disp(T_exact2);

disp('alpha3 ='); disp(alpha3);
disp('T_exact3 ='); disp(T_exact3);

% disp('alpha4 ='); disp(alpha4);
% disp('T_exact4 ='); disp(T_exact4);

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

function Tq = interp_from_aux(aux, tq)
    % Obtiene T(tq) usando 30 puntos finales (excluyendo la última celda rellena)
    % con interpolación polinómica de grado 4.

    % --- Obtener la cell 'save_vars' de forma robusta ---
    if isfield(aux, 'save_vars')
        sv = aux.save_vars;
    elseif isfield(aux, 'aux') && isfield(aux.aux, 'save_vars')
        sv = aux.aux.save_vars;
    else
        error('No se encontró la cell "save_vars" en la estructura proporcionada.');
    end

    % --- Helper inline: escalar numérico finito ---
    isScalarNum = @(x) ~isempty(x) && isnumeric(x) && isscalar(x) && isfinite(x);

    % --- 1) Localizar última columna con valor en fila 2 (temperatura) ---
    lastFilled = [];
    for j = size(sv,2):-1:1
        if isScalarNum(sv{2,j})
            lastFilled = j;
            break
        end
    end
    if isempty(lastFilled)
        error('No hay valores numéricos en la fila 2 (temperatura).');
    end

    % --- 2) Recolectar 5 pares válidos anteriores (tiempo-fila1, temp-fila2) ---
    t_list = []; T_list = [];
    j = lastFilled - 1; % excluir la última rellena
    while numel(T_list) < 50 && j >= 1
        if isScalarNum(sv{1,j}) && isScalarNum(sv{2,j})
            t_list(end+1,1) = sv{1,j}; %#ok<AGROW>
            T_list(end+1,1) = sv{2,j}; %#ok<AGROW>
        end
        j = j - 1;
    end
    if numel(T_list) < 50
        error('No se han encontrado 30 pares válidos anteriores a la última columna rellena.');
    end

    % Orden cronológico y por si acaso ordenar por tiempo ascendente
    t5 = flipud(t_list);
    T5 = flipud(T_list);
    [ts, ord] = sort(t5, 'ascend');
    Ts        = T5(ord);

    % Aviso si tq cae fuera del rango de los 5 puntos (extrapolación)
    if tq < min(ts) || tq > max(ts)
        warning('tq=%.8f cae fuera del rango [%.8f, %.8f]; se extrapola.', tq, min(ts), max(ts));
    end

    % --- 3) Interpolación muy precisa con polinomio de grado 4 ---
    p  = Polyfit(ts, Ts, 6);
    Tq = Polyval(p, tq);
end