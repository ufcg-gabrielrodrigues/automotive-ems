%% Inicializa modelo no Simulink

open_system('models/HybridMPPT.slx', 'loadonly');

%% Parâmetros temporais

T_s = 1.0e-6;   % Passo de cálculo utilizado pelo 'solver' [s]
T_r = 1.0e-5;   % Passo de amostragem global de leitura de variáveis [s]
T_k = 1.0e-4;   % Passo de amostragem global de rotinas de controle [s]
t_f = 5.0e-1;   % Tempo total de simulação [s]

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
i_f_max = 4.5;                  % [A]

% Temperatura da resistência do circuito de estator
alternator.stator.r.T = 150;    % [oC]

%% Retificador

% Parâmetros de controle
rectifier.control.pwm.f_s = 1/T_k;  % Frequência de chaveamento dos PWMs de controle do circuito retificador [Hz]

%% Bateria

battery.v_nom = 80.0;   % [V]

%% Varredura de parâmetros

% Lista de parâmetros a serem varridos individualmente
dynamic_v_l_list = [false true]';   % Atualização dinâmica da tensão de saída para lei de controle
n_r_list = [2000 7500]';        % Velocidade do alternador [rpm]
r_l_list = 1.0;	% Resistência de carga [Ohm]

% Formação das casos de varredura
param_sweep = [];

for index_dynamic_v_l = 1:length(dynamic_v_l_list)
    
    n_r_sweep = [];
    [r_l_dim, ~] = size(r_l_list);
    
    for index_n_r = 1:length(n_r_list)
        n_r_tmp = n_r_list(index_n_r) * ones(r_l_dim, 1);
        n_r_tmp = [n_r_tmp, r_l_list];
        n_r_sweep = [n_r_sweep; n_r_tmp];
    end
    
    [n_r_sweep_dim, ~] = size(n_r_sweep);
    
    dynamic_v_l_tmp = index_dynamic_v_l * ones(n_r_sweep_dim, 1);
    dynamic_v_l_tmp = [dynamic_v_l_tmp, n_r_sweep];
    param_sweep = [param_sweep; dynamic_v_l_tmp];
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
    dynamic_v_l = test_case(1);
    n_r = test_case(2);
    r_l = test_case(3);
    
    simIn(test_case_index) = Simulink.SimulationInput('HybridMPPT');
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('HybridMPPT/dynamic_v_l','Value', num2str(dynamic_v_l));
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
    dynamic_v_l = test_case(1);
    n_r = test_case(2);
    r_l = test_case(3);
    
    % Esquema de controle
    test_case_out(test_case_index).dynamic_v_l = dynamic_v_l;
    
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

%% Carregamento de resultados de an�lise de pot�ncia

try
    simEnv = load('results/PowerAnalysis/simEnv.mat', 'i_f_list', 'n_r_list');
    P_v_o = load('results/PowerAnalysis/P_v_o.mat', 'P_v_o_mpp_sim');
    power_analysis_flag = true;
catch
    power_analysis_flag = false;
end

%% Resultados comparativos

figure_index = 0;

test_cases = param_sweep(param_sweep(:, 1) == 1, 2:3);

for test_case_index = 1:num_cases/length(dynamic_v_l_list)
    % 
    i_f = i_f_max;
    n_r = test_cases(test_case_index, 1);
    r_l = test_cases(test_case_index, 2);
    
    % 
    i_f_index = find(simEnv.i_f_list == i_f);
    n_r_index = find(simEnv.n_r_list == n_r);
    
    % 
    figure_index = figure_index + 1;
	figure(figure_index)
    
    subplot(3, 1, 1)
    
    plot(test_case_out(num_cases/length(dynamic_v_l_list)*0 + test_case_index).electrical_load.p, 'r-');
    hold on;
    plot(test_case_out(num_cases/length(dynamic_v_l_list)*1 + test_case_index).electrical_load.p, 'g-');
    
    if (power_analysis_flag)
        hold on;
        plot([0 t_f], P_v_o.P_v_o_mpp_sim(n_r_index, i_f_index)*[1 1], 'b--');
    end
    
    ylim([0 Inf]);
    title('Pot{\^{e}}ncia el{\''{e}}trica fornecida para a carga');
    xlabel('$t$ $[s]$');
    ylabel('$P_{l}$ $[W]$');
    
    if (power_analysis_flag)
        legend('Tens{\~{a}}o est{\''{a}}tica na carga', 'Atualiza{\c{c}}{\~{a}}o din{\^{a}}mica de tens{\~{a}}o na carga', ...
            'M{\''{a}}xima pot{\^{e}ncia}', 'Location', 'SouthEast');
    else
        legend('Tens{\~{a}}o est{\''{a}}tica na carga', 'Atualiza{\c{c}}{\~{a}}o din{\^{a}}mica de tens{\~{a}}o na carga', ...
            'Location', 'SouthEast');
    end
    
    grid on;
    
    subplot(3, 1, 2)
    plot(test_case_out(num_cases/length(dynamic_v_l_list)*0 + test_case_index).electrical_load.v, 'r-');
    hold on;
    plot(test_case_out(num_cases/length(dynamic_v_l_list)*1 + test_case_index).electrical_load.v, 'g-');
    ylim([0 Inf]);
    title('Tens{\~{a}}o sobre a carga');
    xlabel('$t$ $[s]$');
    ylabel('$v_{l}$ $[V]$');
    legend('Tens{\~{a}}o est{\''{a}}tica na carga', 'Atualiza{\c{c}}{\~{a}}o din{\^{a}}mica de tens{\~{a}}o na carga', ...
        'Location', 'SouthEast');
    grid on;
    
    subplot(3, 1, 3)
    plot(test_case_out(num_cases/length(dynamic_v_l_list)*0 + test_case_index).rectifier.control.u, 'r-');
    hold on;
    plot(test_case_out(num_cases/length(dynamic_v_l_list)*1 + test_case_index).rectifier.control.u, 'g-');
    ylim([0 Inf]);
    title('Ciclo de trabalho aplicado ao retificador');
    xlabel('$t$ $[s]$');
    ylabel('$d_{smr}$');
    legend('Tens{\~{a}}o est{\''{a}}tica na carga', 'Atualiza{\c{c}}{\~{a}}o din{\^{a}}mica de tens{\~{a}}o na carga', ...
        'Location', 'SouthEast');
    grid on;
    
    suptitle(['Caso de teste: $n_{r} = ' num2str(n_r) '$ $[rpm]$; $r_{l} = ' num2str(r_l) '$ $[\Omega]$']);
end

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = sprintf('results/HybridMPPT/Figura_%d', i);
    saveFigure(figure(i), fileName, 'fig');
end

%% Armazenamento dos resultados de simulação

save('results/HybridMPPT/test_case_out.mat', 'test_case_out', '-v7.3');
