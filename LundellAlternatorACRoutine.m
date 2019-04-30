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
% [velocidade do rotor | corrente de excitação]
param_sweep = [1997 1.25;
               1993 1.50;
               2015 1.75;
               3997 1.00;
               3994 1.25;
               3990 1.50;
               3936 1.75;
               5995 1.00;
               5991 1.25;
               5985 1.50;
               5980 1.62;
               5968 1.75];
           
% Resultados experimentais para comparação
% [tensão terminal de linha | corrente terminal de linha]
experimental_results = [7.30 25.17;
                        8.64 29.77;
                        10.07 34.63;
                        7.69 26.50;
                        9.57 33.00;
                        11.42 39.30;
                        13.20 45.35;
                        8.37 28.90;
                        10.26 35.40;
                        12.28 42.30;
                        13.25 45.65;
                        14.28 49.00];

% Resultados simulados para comparação
% [força contra-eletromotriz induzida | tensão terminal de linha | corrente terminal de linha]
simulated_results = [10.59 7.38 25.10;
                     12.47 8.78 29.91;
                     14.26 10.11 34.50;
                     17.21 7.74 26.37;
                     21.31 9.72 33.09;
                     25.04 11.63 39.95;
                     27.93 13.36 45.30;
                     25.88 8.34 28.42;
                     31.78 10.39 35.35;
                     37.30 12.43 42.22;
                     39.66 13.36 45.43;
                     42.07 14.20 48.98];

%% Alternador

% Efeito térmico na resistência do circuito de estator
alternator.stator.r.T = 32; % [oC]

%% Carga elétrica

r_l = 0.509;    % [Ohm]

%% Inicializa modelo no Simulink

open_system('models/LundellAlternatorAC.slx', 'loadonly');

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
set_param('LundellAlternatorAC/Solver Configuration', 'UseLocalSolver', 'on');
set_param('LundellAlternatorAC/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('LundellAlternatorAC/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('LundellAlternatorAC/Solver Configuration', 'DoFixedCost', 'on');
set_param('LundellAlternatorAC/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
set_param('LundellAlternatorAC', 'StopTime', num2str(t_f));

% Salva mundanças feitas no modelo
save_system('models/LundellAlternatorAC.slx');

%% Configuração dos casos de teste como entrada do modelo no Simulink

% Determinação da quantidade de casos de teste
[num_cases, ~] = size(param_sweep);

% Laço de configuração das entradas do modelo no Simulink
for test_case_index = 1:num_cases
    % Isolamento do caso de teste de acordo com a indexação
    test_case = param_sweep(test_case_index, :);
    n_r = test_case(1);
    i_f = test_case(2);
    
    % Configuração da entradas do modelo no Simulink de acordo com a indexação
    simIn(test_case_index) = Simulink.SimulationInput('LundellAlternatorAC');
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('LundellAlternatorAC/n_r', 'Value', num2str(n_r));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('LundellAlternatorAC/i_f', 'Value', num2str(i_f));
end

%% Execução da simulação em ambiente Simulink

% Execução da simulação paralelizada
simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
    'TransferBaseWorkspaceVariables', 'on');

%% Salva e finaliza modelo no Simulink

save_system('models/LundellAlternatorAC.slx');
close_system('models/LundellAlternatorAC.slx');

%% Registro de resultados obtidos na simulação

% Laço de iterações por casos de teste
for test_case_index = 1:num_cases
    % Caso de teste
    test_case_out(test_case_index).test_case = param_sweep(test_case_index, :);
    n_r = test_case_out(test_case_index).test_case(1);
    i_f = test_case_out(test_case_index).test_case(2);
    
    % Alternador
    test_case_out(test_case_index).alternator.rotor.n = n_r;
    test_case_out(test_case_index).alternator.rotor.l.i = simOut(test_case_index).i_f;

    test_case_out(test_case_index).alternator.stator.input.e = simOut(test_case_index).e_a_abc;
    test_case_out(test_case_index).alternator.stator.output.v = simOut(test_case_index).v_a_abc;
    test_case_out(test_case_index).alternator.stator.output.i = simOut(test_case_index).i_a_abc;
    
    test_case_out(test_case_index).alternator.stator.output.e_ll = simOut(test_case_index).e_ll_rms;
    test_case_out(test_case_index).alternator.stator.output.v_ll = simOut(test_case_index).v_ll_rms;
    test_case_out(test_case_index).alternator.stator.output.i_l = simOut(test_case_index).i_l_rms;

    % Carga
    test_case_out(test_case_index).electrical_load.r = r_l;
end

%% Traço de resultados

% Inicialização do índice de figuras
figure_index = 0;

% Laço de iterações por casos de teste
for test_case_index = 1:num_cases
    
    figure_index = figure_index + 1;
    
    figure(figure_index)
    subplot(3, 1, 1)
    plot(test_case_out(test_case_index).alternator.stator.output.e_ll.time, test_case_out(test_case_index).alternator.stator.output.e_ll.data, 'b-');
    hold on;
    plot([0 t_f], [simulated_results(test_case_index, 1) simulated_results(test_case_index, 1)], 'g--');
    hold on;
    error_sim = abs(simulated_results(test_case_index, 1) - test_case_out(test_case_index).alternator.stator.output.e_ll.data);
    plot(test_case_out(test_case_index).alternator.stator.output.e_ll.time, error_sim, 'g:');
    ylabel('$e_{ll}^{RMS} [V]$');
    ylim([0 1.05*max([test_case_out(test_case_index).alternator.stator.output.e_ll.data(end) simulated_results(test_case_index, 1)])]);
    legend('Resultado simulado', 'Resultado simulado referenciado', ...
        'Erro absoluto relativo ao resultado simulado referenciado', 'Location', 'SouthEast');
    grid on;
    
    subplot(3, 1, 2)
    plot(test_case_out(test_case_index).alternator.stator.output.v_ll.time, test_case_out(test_case_index).alternator.stator.output.v_ll.data, 'b-');
    hold on;
    plot([0 t_f], [experimental_results(test_case_index, 1) experimental_results(test_case_index, 1)], 'r--');
    hold on;
    plot([0 t_f], [simulated_results(test_case_index, 2) simulated_results(test_case_index, 2)], 'g--');
    hold on;
    error_exp = abs(experimental_results(test_case_index, 1) - test_case_out(test_case_index).alternator.stator.output.v_ll.data);
    plot(test_case_out(test_case_index).alternator.stator.output.v_ll.time, error_exp, 'r:');
    hold on;
    error_sim = abs(simulated_results(test_case_index, 2) - test_case_out(test_case_index).alternator.stator.output.v_ll.data);
    plot(test_case_out(test_case_index).alternator.stator.output.v_ll.time, error_sim, 'g:');
    ylabel('$v_{ll}^{RMS} [V]$');
    ylim([0 1.05*max([test_case_out(test_case_index).alternator.stator.output.v_ll.data(end) experimental_results(test_case_index, 1) simulated_results(test_case_index, 2)])]);
    legend('Resultado simulado', 'Resultado experimental referenciado', ...
        'Resultado simulado referenciado', ...
        'Erro absoluto relativo ao resultado experimental referenciado', ...
        'Erro absoluto relativo ao resultado simulado referenciado', 'Location', 'SouthEast');
    grid on;
    
    subplot(3, 1, 3)
    plot(test_case_out(test_case_index).alternator.stator.output.i_l.time, test_case_out(test_case_index).alternator.stator.output.i_l.data, 'b-');
    hold on;
    plot([0 t_f], [experimental_results(test_case_index, 2) experimental_results(test_case_index, 2)], 'r--');
    hold on;
    plot([0 t_f], [simulated_results(test_case_index, 3) simulated_results(test_case_index, 3)], 'g--');
    hold on;
    error_exp = abs(experimental_results(test_case_index, 2) - test_case_out(test_case_index).alternator.stator.output.i_l.data);
    plot(test_case_out(test_case_index).alternator.stator.output.i_l.time, error_exp, 'r:');
    hold on;
    error_sim = abs(simulated_results(test_case_index, 3) - test_case_out(test_case_index).alternator.stator.output.i_l.data);
    plot(test_case_out(test_case_index).alternator.stator.output.i_l.time, error_sim, 'g:');
    xlabel('$t [s]$');
    ylabel('$i_{l}^{RMS} [A]$');
    ylim([0 1.05*max([test_case_out(test_case_index).alternator.stator.output.i_l.data(end) experimental_results(test_case_index, 2) simulated_results(test_case_index, 3)])]);
    legend('Resultado simulado', 'Resultado experimental referenciado', ...
        'Resultado simulado referenciado', ...
        'Erro absoluto relativo ao resultado experimental referenciado', ...
        'Erro absoluto relativo ao resultado simulado referenciado', 'Location', 'SouthEast');
    grid on;
end

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = sprintf('results/LundellAlternator/alternador-validacao-ca-%d', i);
    saveFigure(figure(i), fileName, 'fig');
    saveFigure(figure(i), fileName, 'png');
end

%% 

% 
e_ll_validation = zeros(num_cases, 5);
v_ll_validation = zeros(num_cases, 5);
i_l_validation = zeros(num_cases, 5);

% 
for test_case_index = 1:num_cases
    e_ll_validation(test_case_index, 1) = mean(test_case_out(test_case_index).alternator.stator.output.e_ll.data(ceil(end/2), end));
    e_ll_validation(test_case_index, 2) = NaN;
    e_ll_validation(test_case_index, 3) = NaN;
    e_ll_validation(test_case_index, 4) = simulated_results(test_case_index, 1);
    e_ll_validation(test_case_index, 5) = 100*abs(e_ll_validation(test_case_index, 1) - e_ll_validation(test_case_index, 4))/e_ll_validation(test_case_index, 4);
    
    v_ll_validation(test_case_index, 1) = mean(test_case_out(test_case_index).alternator.stator.output.v_ll.data(ceil(end/2), end));
    v_ll_validation(test_case_index, 2) = experimental_results(test_case_index, 1);
    v_ll_validation(test_case_index, 3) = 100*abs(v_ll_validation(test_case_index, 1) - v_ll_validation(test_case_index, 2))/v_ll_validation(test_case_index, 2);
    v_ll_validation(test_case_index, 4) = simulated_results(test_case_index, 2);
    v_ll_validation(test_case_index, 5) = 100*abs(v_ll_validation(test_case_index, 1) - v_ll_validation(test_case_index, 4))/v_ll_validation(test_case_index, 4);
    
    i_l_validation(test_case_index, 1) = mean(test_case_out(test_case_index).alternator.stator.output.i_l.data(ceil(end/2), end));
    i_l_validation(test_case_index, 2) = experimental_results(test_case_index, 2);
    i_l_validation(test_case_index, 3) = 100*abs(i_l_validation(test_case_index, 1) - i_l_validation(test_case_index, 2))/i_l_validation(test_case_index, 2);
    i_l_validation(test_case_index, 4) = simulated_results(test_case_index, 3);
    i_l_validation(test_case_index, 5) = 100*abs(i_l_validation(test_case_index, 1) - i_l_validation(test_case_index, 4))/i_l_validation(test_case_index, 4);
end

% 
xlswrite('results/LundellAlternator/validacao_ca.xls', e_ll_validation, 'e_ll');
xlswrite('results/LundellAlternator/validacao_ca.xlsx', v_ll_validation, 'v_ll');
xlswrite('results/LundellAlternator/validacao_ca.xlsx', i_l_validation, 'i_l');

%% Armazenamento dos resultados de simulação

save('results/LundellAlternator/test_case_out_ac.mat', 'test_case_out', '-v7.3');

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
