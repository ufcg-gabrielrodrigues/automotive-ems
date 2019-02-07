%% Inicializa modelo no Simulink

open_system('models/HybridMPPT.slx', 'loadonly');

%% Parâmetros temporais

T_s = 1.0e-6;   % Passo de cálculo utilizado pelo 'solver' [s]
T_k = 1.0e-4;   % Passo de amostragem global de rotinas de controle [s]
t_f = 1.0e-0;   % Tempo total de simulação [s]

%% Esquema de controle

% Atualização de parâmetro: fator de acoplamento
k_e_default_local = 'k_e = 0;';

if (isfield(alternator.k_e, 'function'))
    k_e_str = regexprep(func2str(alternator.k_e.function), '@\(.+?\)', '');
else
    k_e_str = num2str(alternator.k_e.value);
end

if (alternator.stator.connection == delta)
    k_e_str = ['(' k_e_str ')./sqrt(3)'];
end

k_e_local = ['k_e = ' k_e_str ';'];
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'HybridMPPT/Load Matching Control [Voltage-based]/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, k_e_default_local, k_e_local);

%% Alternador

% Corrente de excitação máxima
i_f_max = 4.5;  % [A]

% Efeito térmico na resistência do circuito de estator
T = 150;        % [oC]
alternator.stator.r.value = alternator.stator.r.function(T);

%% Bateria

battery.v_nom = 80.0;   % [V]

%% Varredura de parâmetros

% Lista de parâmetros a serem varridos individualmente
dynamic_v_o_list = [false true]';   % Atualização dinâmica da tensão de saída para lei de controle
n_r_list = (2000:500:7500)';        % Velocidade do alternador [rpm]
r_l_list = [0.15 (0.5:0.5:2.0)]';	% Resistência de carga [Ohm]

% Formação das casos de varredura
param_sweep = [];

for index_dynamic_v_o = 1:length(dynamic_v_o_list)
    
    n_r_sweep = [];
    [r_l_dim, ~] = size(r_l_list);
    
    for index_n_r = 1:length(n_r_list)
        n_r_tmp = n_r_list(index_n_r) * ones(r_l_dim, 1);
        n_r_tmp = [n_r_tmp, r_l_list];
        n_r_sweep = [n_r_sweep; n_r_tmp];
    end
    
    [n_r_sweep_dim, ~] = size(n_r_sweep);
    
    dynamic_v_o_tmp = index_dynamic_v_o * ones(n_r_sweep_dim, 1);
    dynamic_v_o_tmp = [dynamic_v_o_tmp, n_r_sweep];
    param_sweep = [param_sweep; dynamic_v_o_tmp];
end

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
set_param('HybridMPPT/Solver Configuration', 'UseLocalSolver', 'on');
set_param('HybridMPPT/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('HybridMPPT/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('HybridMPPT/Solver Configuration', 'DoFixedCost', 'on');
set_param('HybridMPPT/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
set_param('HybridMPPT', 'StopTime', num2str(t_f));

% Salva mundanças feitas no modelo
save_system('models/HybridMPPT.slx');

%% Configuração dos casos de teste como entrada do modelo no Simulink

[num_cases, ~] = size(param_sweep);

for test_case_index = 1:num_cases
    test_case = param_sweep(test_case_index, :);
    dynamic_v_o = test_case(1);
    n_r = test_case(2);
    r_l = test_case(3);
    
    simIn(test_case_index) = Simulink.SimulationInput('HybridMPPT');
    simIn(test_case_index) = simIn(test_case_index).setVariable('HybridMPPT/dynamic_v_o','Value', num2str(dynamic_v_o));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('HybridMPPT/n_r', 'Value', num2str(n_r));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('HybridMPPT/Load/r_l', 'R', num2str(r_l));
end

%% Execução da simulação em ambiente Simulink

% Execução da simulação paralelizada
simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
    'TransferBaseWorkspaceVariables', 'on');

%% Redefinição de parâmetros de alternador

% Atualização de parâmetro para valor padrão: fator de acoplamento
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'HybridMPPT/Load Matching Control [Voltage-based]/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, k_e_local, k_e_default_local);

%% Finaliza modelo no Simulink

save_system('models/HybridMPPT.slx');
close_system('models/HybridMPPT.slx');

%% Registro de resultados obtidos no caso de teste

% Laço de iterações por casos de teste
for test_case_index = 1:num_cases
    
    test_case = param_sweep(test_case_index, :);
    
    % Valor de variáveis correspondentes ao caso de teste
    dynamic_v_o = test_case(1);
    n_r = test_case(2);
    r_l = test_case(3);
    
    % Esquema de controle
    test_case_out(test_case_index).dynamic_v_o = dynamic_v_o;
    
    % Alternador
    test_case_out(test_case_index).alternator.rotor.n = n_r;
    test_case_out(test_case_index).alternator.rotor.l.i = simOut(test_case_index).i_f;
    
    % Retificador
    test_case_out(test_case_index).rectifier.control.u = simOut(test_case_index).u_smr;
    
    % Carga
    test_case_out(test_case_index).electrical_load.r = r_l;
    test_case_out(test_case_index).electrical_load.v = simOut(test_case_index).v_l;
    test_case_out(test_case_index).electrical_load.i = simOut(test_case_index).i_l;
    test_case_out(test_case_index).electrical_load.p = simOut(test_case_index).p_l;
end

%% Resultados comparativos

figure_index = 0;

test_cases = param_sweep(param_sweep(:, 1) == 1, 2:3);

for test_case_index = 1:num_cases/length(dynamic_v_o_list)
    
    figure_index = figure_index + 1;
	figure(figure_index)
    
    subplot(3, 1, 1)
    
    plot(test_case_out(num_cases/length(dynamic_v_o_list)*0 + test_case_index).electrical_load.p, 'r-');
    hold on;
    plot(test_case_out(num_cases/length(dynamic_v_o_list)*1 + test_case_index).electrical_load.p, 'g-');
    ylim([0 Inf]);
    title('Pot{\^{e}}ncia el{\''{e}}trica fornecida para a carga');
    xlabel('$t$ $[s]$');
    ylabel('$P_{l}$ $[W]$');
    legend('Tens{\~{a}}o est{\''{a}}tica na carga', 'Atualiza{\c{c}}{\~{a}}o din{\^{a}}mica de tens�o na carga', ...
        'Location', 'SouthEast');
    grid on;
    
    subplot(3, 1, 2)
    plot(test_case_out(num_cases/length(dynamic_v_o_list)*0 + test_case_index).electrical_load.v, 'r-');
    hold on;
    plot(test_case_out(num_cases/length(dynamic_v_o_list)*1 + test_case_index).electrical_load.v, 'g-');
    ylim([0 Inf]);
    title('Tens{\~{a}}o sobre a carga');
    xlabel('$t$ $[s]$');
    ylabel('$v_{l}$ $[V]$');
    legend('Tens{\~{a}}o est{\''{a}}tica na carga', 'Atualiza{\c{c}}{\~{a}}o din{\^{a}}mica de tens�o na carga', ...
        'Location', 'SouthEast');
    grid on;
    
    subplot(3, 1, 3)
    plot(test_case_out(num_cases/length(dynamic_v_o_list)*0 + test_case_index).rectifier.control.u, 'r-');
    hold on;
    plot(test_case_out(num_cases/length(dynamic_v_o_list)*1 + test_case_index).rectifier.control.u, 'g-');
    ylim([0 Inf]);
    title('Ciclo de trabalho aplicado ao retificador');
    xlabel('$t$ $[s]$');
    ylabel('$d_{smr}$');
    legend('Tens{\~{a}}o est{\''{a}}tica na carga', 'Atualiza{\c{c}}{\~{a}}o din{\^{a}}mica de tens�o na carga', ...
        'Location', 'SouthEast');
    grid on;
    
    suptitle(['Caso de teste: $n_{r} = ' num2str(test_cases(test_case_index, 1)) ...
        '$ $[rpm]$; $r_{l} = ' num2str(test_cases(test_case_index, 2)) '$ $[\Omega]$']);
end

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = sprintf('results/HybridMPPT/Figura_%d', i);
    saveFigure(figure(i), fileName, 'fig');
end

%% Armazenamento dos resultados de simulação

save('results/HybridMPPT/test_case_out.mat', 'test_case_out', '-v7.3');
