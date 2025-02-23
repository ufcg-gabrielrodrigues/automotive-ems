%% Registro de variáveis existente no workspace antes da execução do script

% Registra variáveis existentes no workspace neste instante
if (exist('previousVars', 'var'))
    previousVars{:, end + 1} = who;
else
    previousVars{1} = who;
end

%% 

close all;

%% Parâmetros temporais

T_s = 1e-6; % Passo de cálculo utilizado pelo 'solver' [s]
t_f = 2e-1; % Tempo total de simulação [s]

%% Varredura de parâmetros

% Parâmetros afixados
% [velocidade do rotor | corrente na carga | tensão na carga]
param_sweep = [1967 36.058 13.49;
               3994 36.240 13.54;
               3983 55.650 13.58;
               5992 36.220 13.61;
               5985 56.190 13.46;
               5967 75.610 13.58];
           
% Resultados experimentais para comparação
% [corrente de excitação | tensão terminal de linha | corrente terminal de linha]
experimental_results = [1.970 12.26 27.90;
                        1.250 12.42 27.50;
                        1.650 12.55 41.70;
                        1.060 12.53 27.30;
                        1.510 12.46 42.00;
                        1.970 12.60 56.30];

% Resultados simulados para comparação
% [corrente de excitação | tensão terminal de linha | corrente terminal de linha]
simulated_results = [1.960 12.70 27.20;
                     1.199 12.75 26.90;
                     1.610 12.95 41.30;
                     1.034 12.77 26.87;
                     1.490 12.80 41.60;
                     1.950 13.00 56.00];

%% Alternador

% Referência de tensão CC
v_dc_ref = 13.5;            % [V]

% Corrente de excitação máxima
i_f_max = 5.0e-0;           % [A]

% Efeito térmico na resistência do circuito de estator
alternator.stator.r.T = 32; % [oC]

%% Retificador

% Filtro passivo
rectifier.filter.c = 1e-3;  % Capacitância de filtro [F]

%% Inicializa modelo no Simulink

open_system('models/LundellAlternatorDC.slx', 'loadonly');

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
set_param('LundellAlternatorDC/Solver Configuration', 'UseLocalSolver', 'on');
set_param('LundellAlternatorDC/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('LundellAlternatorDC/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('LundellAlternatorDC/Solver Configuration', 'DoFixedCost', 'on');
set_param('LundellAlternatorDC/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
set_param('LundellAlternatorDC', 'StopTime', num2str(t_f));

% Salva mundanças feitas no modelo
save_system('models/LundellAlternatorDC.slx');

%% Configuração dos casos de teste como entrada do modelo no Simulink

% Determinação da quantidade de casos de teste
[num_cases, ~] = size(param_sweep);

% Laço de configuração das entradas do modelo no Simulink
for test_case_index = 1:num_cases
    % Isolamento do caso de teste de acordo com a indexação
    test_case = param_sweep(test_case_index, :);
    n_r = test_case(1);
    i_dc = test_case(2);
    v_dc = test_case(3);
    
    % Determinação da resistência de carga
    r_dc = v_dc/i_dc;
    
    % Configuração da entradas do modelo no Simulink de acordo com a indexação
    simIn(test_case_index) = Simulink.SimulationInput('LundellAlternatorDC');
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('LundellAlternatorDC/n_r', 'Value', num2str(n_r));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('LundellAlternatorDC/r_l', 'R', num2str(r_dc));
end

%% Execução da simulação em ambiente Simulink

% Execução da simulação paralelizada
simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
    'TransferBaseWorkspaceVariables', 'on');

%% Salva e finaliza modelo no Simulink

save_system('models/LundellAlternatorDC.slx');
close_system('models/LundellAlternatorDC.slx');

%% Registro de resultados obtidos na simulação

% Laço de iterações por casos de teste
for test_case_index = 1:num_cases
    % Caso de teste
    test_case_out(test_case_index).test_case = param_sweep(test_case_index, :);
    n_r = test_case_out(test_case_index).test_case(1);
    i_dc = test_case_out(test_case_index).test_case(2);
    v_dc = test_case_out(test_case_index).test_case(3);
    
    % Alternador
    test_case_out(test_case_index).alternator.rotor.n = n_r;
    test_case_out(test_case_index).alternator.rotor.l.i = simOut(test_case_index).i_f;

    test_case_out(test_case_index).alternator.stator.input.e = simOut(test_case_index).e_a_abc;
    test_case_out(test_case_index).alternator.stator.output.v = simOut(test_case_index).v_a_abc;
    test_case_out(test_case_index).alternator.stator.output.i = simOut(test_case_index).i_a_abc;
    
    test_case_out(test_case_index).alternator.stator.output.e_ll = simOut(test_case_index).e_ll_rms;
    test_case_out(test_case_index).alternator.stator.output.v_ll = simOut(test_case_index).v_ll_rms;
    test_case_out(test_case_index).alternator.stator.output.i_l = simOut(test_case_index).i_l_rms;

    % Bateria
    test_case_out(test_case_index).battery.v = simOut(test_case_index).v_b;

    % Carga
    test_case_out(test_case_index).electrical_load.r = v_dc/i_dc;
    test_case_out(test_case_index).electrical_load.v = simOut(test_case_index).v_l;
    test_case_out(test_case_index).electrical_load.i = simOut(test_case_index).i_l;
    test_case_out(test_case_index).electrical_load.p = simOut(test_case_index).p_l;
end

%% Traço de resultados

if (plotDCValidationFlag)
    % Inicialização do índice de figuras
    figure_index = 0;
    
    % Cópia do vetor de tempo da simulação
    time = test_case_out(test_case_index).alternator.rotor.l.i.time;
    
    % Laço de iterações por casos de teste
    for test_case_index = 1:num_cases
        
        figure_index = figure_index + 1;
        
        figure(figure_index)
        subplot(3, 1, 1)
        plot(test_case_out(test_case_index).alternator.rotor.l.i.time, test_case_out(test_case_index).alternator.rotor.l.i.data, 'b-');
        hold on;
        plot([0 t_f], [experimental_results(test_case_index, 1) experimental_results(test_case_index, 1)], 'r--');
        hold on;
        plot([0 t_f], [simulated_results(test_case_index, 1) simulated_results(test_case_index, 1)], 'g--');
        hold on;
        error_exp = abs(experimental_results(test_case_index, 1) - test_case_out(test_case_index).alternator.rotor.l.i.data);
        plot(time, error_exp, 'r:');
        hold on;
        error_sim = abs(simulated_results(test_case_index, 1) - test_case_out(test_case_index).alternator.rotor.l.i.data);
        plot(time, error_sim, 'g:');
        ylabel('$i_{f}\,[\textrm{A}]$');
        ylim([0 1.05*max([test_case_out(test_case_index).alternator.rotor.l.i.data(end) experimental_results(test_case_index, 1) simulated_results(test_case_index, 1)])]);
        legend('Resultado simulado', 'Resultado experimental referenciado', ...
            'Resultado simulado referenciado', ...
            'Erro absoluto relativo ao resultado experimental referenciado', ...
            'Erro absoluto relativo ao resultado simulado referenciado', 'Location', 'SouthEast');
        grid on;
        
        subplot(3, 1, 2)
        plot(test_case_out(test_case_index).alternator.stator.output.v_ll.time, test_case_out(test_case_index).alternator.stator.output.v_ll.data, 'b-');
        hold on;
        plot([0 t_f], [experimental_results(test_case_index, 2) experimental_results(test_case_index, 2)], 'r--');
        hold on;
        plot([0 t_f], [simulated_results(test_case_index, 2) simulated_results(test_case_index, 2)], 'g--');
        hold on;
        error_exp = abs(experimental_results(test_case_index, 2) - test_case_out(test_case_index).alternator.stator.output.v_ll.data);
        plot(time, error_exp, 'r:');
        hold on;
        error_sim = abs(simulated_results(test_case_index, 2) - test_case_out(test_case_index).alternator.stator.output.v_ll.data);
        plot(time, error_sim, 'g:');
        ylabel('$v_{s,\textrm{rms}}^{ll}\,[\textrm{V}]$');
        ylim([0 1.05*max([test_case_out(test_case_index).alternator.stator.output.v_ll.data(end) experimental_results(test_case_index, 2) simulated_results(test_case_index, 2)])]);
        legend('Resultado simulado', 'Resultado experimental referenciado', ...
            'Resultado simulado referenciado', ...
            'Erro absoluto relativo ao resultado experimental referenciado', ...
            'Erro absoluto relativo ao resultado simulado referenciado', 'Location', 'SouthEast');
        grid on;
        
        subplot(3, 1, 3)
        plot(test_case_out(test_case_index).alternator.stator.output.i_l.time, test_case_out(test_case_index).alternator.stator.output.i_l.data, 'b-');
        hold on;
        plot([0 t_f], [experimental_results(test_case_index, 3) experimental_results(test_case_index, 3)], 'r--');
        hold on;
        plot([0 t_f], [simulated_results(test_case_index, 3) simulated_results(test_case_index, 3)], 'g--');
        hold on;
        error_exp = abs(experimental_results(test_case_index, 3) - test_case_out(test_case_index).alternator.stator.output.i_l.data);
        plot(time, error_exp, 'r:');
        hold on;
        error_sim = abs(simulated_results(test_case_index, 3) - test_case_out(test_case_index).alternator.stator.output.i_l.data);
        plot(time, error_sim, 'g:');
        xlabel('$t [s]$');
        ylabel('$i_{s,\textrm{rms}}^{ll}\,[\textrm{A}]$');
        ylim([0 1.05*max([test_case_out(test_case_index).alternator.stator.output.i_l.data(end) experimental_results(test_case_index, 3) simulated_results(test_case_index, 3)])]);
        legend('Resultado simulado', 'Resultado experimental referenciado', ...
            'Resultado simulado referenciado', ...
            'Erro absoluto relativo ao resultado experimental referenciado', ...
            'Erro absoluto relativo ao resultado simulado referenciado', 'Location', 'SouthEast');
        grid on;
    end
    
    % Armazenamento de figuras
    for i = 1:figure_index
        fileName = sprintf('results/LundellAlternator/alternador-validacao-cc-%d', i);
        saveFigure(figure(i), fileName, 'fig');
        saveFigure(figure(i), fileName, 'png');
    end
end

%% 

% 
i_f_validation = zeros(num_cases, 5);
v_ll_validation = zeros(num_cases, 5);
i_l_validation = zeros(num_cases, 5);

% 
for test_case_index = 1:num_cases
    i_f_validation(test_case_index, 1) = mean(test_case_out(test_case_index).alternator.rotor.l.i.data(ceil(end/2), end));
    i_f_validation(test_case_index, 2) = experimental_results(test_case_index, 1);
    i_f_validation(test_case_index, 3) = 100*abs(i_f_validation(test_case_index, 1) - i_f_validation(test_case_index, 2))/i_f_validation(test_case_index, 2);
    i_f_validation(test_case_index, 4) = simulated_results(test_case_index, 1);
    i_f_validation(test_case_index, 5) = 100*abs(i_f_validation(test_case_index, 1) - i_f_validation(test_case_index, 4))/i_f_validation(test_case_index, 4);
    
    v_ll_validation(test_case_index, 1) = mean(test_case_out(test_case_index).alternator.stator.output.v_ll.data(ceil(end/2), end));
    v_ll_validation(test_case_index, 2) = experimental_results(test_case_index, 2);
    v_ll_validation(test_case_index, 3) = 100*abs(v_ll_validation(test_case_index, 1) - v_ll_validation(test_case_index, 2))/v_ll_validation(test_case_index, 2);
    v_ll_validation(test_case_index, 4) = simulated_results(test_case_index, 2);
    v_ll_validation(test_case_index, 5) = 100*abs(v_ll_validation(test_case_index, 1) - v_ll_validation(test_case_index, 4))/v_ll_validation(test_case_index, 4);
    
    i_l_validation(test_case_index, 1) = mean(test_case_out(test_case_index).alternator.stator.output.i_l.data(ceil(end/2), end));
    i_l_validation(test_case_index, 2) = experimental_results(test_case_index, 3);
    i_l_validation(test_case_index, 3) = 100*abs(i_l_validation(test_case_index, 1) - i_l_validation(test_case_index, 2))/i_l_validation(test_case_index, 2);
    i_l_validation(test_case_index, 4) = simulated_results(test_case_index, 3);
    i_l_validation(test_case_index, 5) = 100*abs(i_l_validation(test_case_index, 1) - i_l_validation(test_case_index, 4))/i_l_validation(test_case_index, 4);
end

% 
xlswrite('results/LundellAlternator/validacao_cc.xlsx', [param_sweep i_f_validation], 'i_f');
xlswrite('results/LundellAlternator/validacao_cc.xlsx', [param_sweep v_ll_validation], 'v_ll');
xlswrite('results/LundellAlternator/validacao_cc.xlsx', [param_sweep i_l_validation], 'i_l');

%% Armazenamento dos resultados de simulação

save('results/LundellAlternator/test_case_out_dc.mat', 'test_case_out', '-v7.3');

%% Exclusão das variáveis excedentes

% Inicializa variáveis para que apareçam no workspace
currentVars = [];
newVars = [];

% Identifica variáveis existentes no workspace neste instante
currentVars = who;

% Pela diferença, determina variáveis criadas pelo script e, se necessário,
% apaga registro atual de variáveis
newVars = setdiff(currentVars, previousVars{end});

% Exclui último registro de variáveis do workspace
if (length(previousVars) > 1)
    previousVars(end) = [];
end

% Realiza limpeza
clear(newVars{:});
