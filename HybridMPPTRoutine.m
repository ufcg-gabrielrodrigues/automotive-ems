%% Inicializa modelo no Simulink

open_system('models/HybridMPPT.slx', 'loadonly');

%% Parâmetros temporais

T_s = 1.0e-6;   % Passo de cálculo utilizado pelo 'solver' [s]
T_k = 1.0e-4;   % Passo de amostragem global de rotinas de controle [s]
t_f = 2.0e-2;   % Tempo total de simulação [s]

%% Esquemas de controle

% Atualização dinâmica da tensão de saída para lei de controle
dynamic_voltage_flag = true;

% Lista de títulos por esquema de controle
% control_scheme_title = {'Esquema de controle baseado em casamento de tens{\~{a}o}', ...
%     'Esquema de controle baseado em casamento de imped{\^{a}}ncia', ...
%     'Esquema de controle baseado em casamento de tens{\~{a}o} ajustado', ...
%     'Esquema de controle baseado em casamento de imped{\^{a}}ncia ajustado'};
control_scheme_title = {'Esquema de controle baseado em casamento de tens{\~{a}o}', ...
    'Esquema de controle baseado em casamento de tens{\~{a}o} ajustado'};

%
voltageMatchingControlScheme = Simulink.Variant('control_scheme == 1');

% 
impedanceMatchingControlScheme = Simulink.Variant('control_scheme == 2');

% 
fittedVoltageMatchingControlScheme = Simulink.Variant('control_scheme == 3');

try
    load('results/PowerAnalysis/P_v_o.mat', 'v_o_mpp_fit');
catch
    disp('Não foi possível carregar a variável v_o_mpp_fit');
end

v_x_default = 'v_x = v_o;';
v_x = ['v_x = ' fitToString(v_o_mpp_fit) ';'];

blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'HybridMPPT/Control scheme/Load Matching Control [Fitted Voltage-based]/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, v_x_default, v_x);

% 
fittedImpedanceMatchingControlScheme = Simulink.Variant('control_scheme == 4');

try
    load('results/PowerAnalysis/P_z_o.mat', 'z_o_mpp_fit');
catch
    disp('Não foi possível carregar a variável z_o_mpp_fit');
end

z_x_default = 'z_x = z_o;';
z_x = ['z_x = ' fitToString(z_o_mpp_fit) ';'];

blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'HybridMPPT/Control scheme/Load Matching Control [Fitted Impedance-based]/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, z_x_default, z_x);

%% Alternador

% Corrente de excitação máxima
i_f_max = 4.5;  % [A]

% Efeito térmico na resistência do circuito de estator
T = 150;        % [oC]
alternator.stator.r.value = alternator.stator.r.function(T);

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
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'HybridMPPT/Control scheme/Load Matching Control [Voltage-based]/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, k_e_default_local, k_e_local);

% Atualização de parâmetro: indutância própria de estator
l_s_default_local = 'l_s = 1e-6;';

if (isfield(alternator.stator.l, 'function'))
    l_s_str = regexprep(func2str(alternator.stator.l.function), '@\(.+?\)', '');
else
    l_s_str = num2str(alternator.stator.l.value);
end

if (alternator.stator.connection == delta)
    l_s_str = ['(' l_s_str ')./3'];
end

l_s_local = ['l_s = ' l_s_str ';'];
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'HybridMPPT/Control scheme/Load Matching Control [Impedance-based]/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, l_s_default_local, l_s_local);

%% Bateria

battery.v_nom = 80.0;   % [V]

%% Varredura de parâmetros

% Lista de parâmetros a serem varridos individualmente
n_r_list = (2000:500:7500)';    % Velocidade do alternador [rpm]
r_l_list = [0.15 (0.5:0.5:2.0)]';      % Resistência de carga [Ohm]

% Formação das casos de varredura
param_sweep = [];

% for index_control_scheme = 1:length(control_scheme_title)
for index_control_scheme = [1 3]
    
    n_r_sweep = [];
    [r_l_dim, ~] = size(r_l_list);
    
    for index_n_r = 1:length(n_r_list)
        n_r_tmp = n_r_list(index_n_r) * ones(r_l_dim, 1);
        n_r_tmp = [n_r_tmp, r_l_list];
        n_r_sweep = [n_r_sweep; n_r_tmp];
    end
    
    [n_r_sweep_dim, ~] = size(n_r_sweep);
    
    control_scheme_tmp = index_control_scheme * ones(n_r_sweep_dim, 1);
    control_scheme_tmp = [control_scheme_tmp, n_r_sweep];
    param_sweep = [param_sweep; control_scheme_tmp];
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

% Inicialização da variável que define o esquema de controle selecionado
control_scheme = 1;

[num_cases, ~] = size(param_sweep);

for test_case_index = 1:num_cases
    test_case = param_sweep(test_case_index, :);
    control_scheme = test_case(1);
    n_r = test_case(2);
    r_l = test_case(3);
    
    simIn(test_case_index) = Simulink.SimulationInput('HybridMPPT');
    simIn(test_case_index) = simIn(test_case_index).setVariable('control_scheme', control_scheme);
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('HybridMPPT/n_r', 'Value', num2str(n_r));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('HybridMPPT/Load/r_l', 'R', num2str(r_l));
end

%% Execução da simulação em ambiente Simulink

% Execução da simulação paralelizada
simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
    'TransferBaseWorkspaceVariables', 'on');

%% Redefinição de parâmetros de alternador

% Atualização de parâmetro para valor padrão: fator de acoplamento
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'HybridMPPT/Control scheme/Load Matching Control [Voltage-based]/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, k_e_local, k_e_default_local);

% Atualização de parâmetro para valor padrão: indutância própria de estator
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'HybridMPPT/Control scheme/Load Matching Control [Impedance-based]/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, l_s_local, l_s_default_local);

% 
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'HybridMPPT/Control scheme/Load Matching Control [Fitted Voltage-based]/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, v_x, v_x_default);

% 
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'HybridMPPT/Control scheme/Load Matching Control [Fitted Impedance-based]/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, z_x, z_x_default);

%% Finaliza modelo no Simulink

save_system('models/HybridMPPT.slx');
close_system('models/HybridMPPT.slx');

%% Registro de resultados obtidos no caso de teste

% Laço de iterações por casos de teste
for test_case_index = 1:num_cases
    
    test_case = param_sweep(test_case_index, :);
    
    % Valor de variáveis correspondentes ao caso de teste
    control_scheme = test_case(1);
    n_r = test_case(2);
    r_l = test_case(3);
    
    % Esquema de controle
    test_case_out(test_case_index).control_scheme = control_scheme;
    
    % Alternador
    test_case_out(test_case_index).alternator.rotor.n = n_r;
    test_case_out(test_case_index).alternator.rotor.l.i = simOut(test_case_index).i_f;
    
%     test_case_out(test_case_index).alternator.stator.input.e.value = simOut(test_case_index).e_a_abc;
%     test_case_out(test_case_index).alternator.stator.output.v = simOut(test_case_index).v_a_abc;
%     test_case_out(test_case_index).alternator.stator.output.i = simOut(test_case_index).i_a_abc;
    
    % Retificador
    test_case_out(test_case_index).rectifier.control.u = simOut(test_case_index).u_smr;
    
    % Carga
    test_case_out(test_case_index).electrical_load.r = r_l;
    test_case_out(test_case_index).electrical_load.v = simOut(test_case_index).v_l;
    test_case_out(test_case_index).electrical_load.i = simOut(test_case_index).i_l;
    test_case_out(test_case_index).electrical_load.z = simOut(test_case_index).z_l;
    test_case_out(test_case_index).electrical_load.p = simOut(test_case_index).p_l;
end

%% Resultados comparativos

figure_index = 0;

test_cases = param_sweep(param_sweep(:, 1) == 1, 2:3);

for test_case_index = 1:num_cases/length(control_scheme_title)
    
    figure_index = figure_index + 1;
	figure(figure_index)
    
    subplot(2, 1, 1)
    
    plot(test_case_out(num_cases/length(control_scheme_title)*0 + test_case_index).electrical_load.p, 'r-');
    hold on;
    plot(test_case_out(num_cases/length(control_scheme_title)*1 + test_case_index).electrical_load.p, 'g-');
%     hold on;
%     plot(test_case_out(num_cases/length(control_scheme_title)*2 + test_case_index).electrical_load.p, 'b-');
%     hold on;
%     plot(test_case_out(num_cases/length(control_scheme_title)*3 + test_case_index).electrical_load.p, 'k-');
    ylim([0 Inf]);
    title(['Pot{\^{e}}ncia el{\''{e}}trica fornecida para a carga ($n_{r} = ' ...
        num2str(test_cases(test_case_index, 1)) '$ $[rpm]$; $r_{l} = ' ...
        num2str(test_cases(test_case_index, 2)) '$ $[\Omega]$)']);
    xlabel('$t$ $[s]$');
    ylabel('$P_{l}$ $[W]$');
%     legend('Casamento de tens{\~{a}}o', 'Casamento de imped{\^{a}}ncia', ...
%         'Casamento de tens{\~{a}}o ajustado', 'Casamento de imped{\^{a}}ncia ajustado', ...
%         'Location', 'NorthEast');
    legend('Casamento de tens{\~{a}}o', 'Casamento de tens{\~{a}}o ajustado', 'Location', 'SouthEast');
    grid on;
    
    subplot(2, 1, 2)
    plot(test_case_out(num_cases/length(control_scheme_title)*0 + test_case_index).electrical_load.v, 'r-');
    hold on;
    plot(test_case_out(num_cases/length(control_scheme_title)*1 + test_case_index).electrical_load.v, 'g-');
%     hold on;
%     plot(test_case_out(num_cases/length(control_scheme_title)*2 + test_case_index).electrical_load.v, 'b-');
%     hold on;
%     plot(test_case_out(num_cases/length(control_scheme_title)*3 + test_case_index).electrical_load.v, 'k-');
    ylim([0 Inf]);
    title(['Tens{\~{a}}o sobre a carga ($n_{r} = ' num2str(test_cases(test_case_index, 1)) ...
        '$ $[rpm]$; $r_{l} = ' num2str(test_cases(test_case_index, 2)) '$ $[\Omega]$)']);
    xlabel('$t$ $[s]$');
    ylabel('$v_{l}$ $[V]$');
%     legend('Casamento de tens{\~{a}}o', 'Casamento de imped{\^{a}}ncia', ...
%         'Casamento de tens{\~{a}}o ajustado', 'Casamento de imped{\^{a}}ncia ajustado', ...
%         'Location', 'NorthEast');
    legend('Casamento de tens{\~{a}}o', 'Casamento de tens{\~{a}}o ajustado', 'Location', 'SouthEast');
    grid on;
    
%     subplot(3, 1, 3)
%     plot(test_case_out(num_cases/length(control_scheme_title)*0 + test_case_index).electrical_load.z, 'r-');
%     hold on;
%     plot(test_case_out(num_cases/length(control_scheme_title)*1 + test_case_index).electrical_load.z, 'g-');
%     hold on;
%     plot(test_case_out(num_cases/length(control_scheme_title)*2 + test_case_index).electrical_load.z, 'b-');
%     hold on;
%     plot(test_case_out(num_cases/length(control_scheme_title)*3 + test_case_index).electrical_load.z, 'k-');
%     title(['Imped{\^{a}}ncia de carga observada ($n_{r} = ' num2str(test_cases(test_case_index, 1)) ...
%         '$ $[rpm]$; $r_{l} = ' num2str(test_cases(test_case_index, 2)) '$ $[\Omega]$)']);
%     xlabel('$t$ $[s]$');
%     ylabel('$z_{l}$ $[\Omega]$');
%     legend('Casamento de tens{\~{a}}o', 'Casamento de imped{\^{a}}ncia', ...
%         'Casamento de tens{\~{a}}o ajustado', 'Casamento de imped{\^{a}}ncia ajustado', ...
%         'Location', 'NorthEast');
%     grid on;
end

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = sprintf('results/HybridMPPT/Figura_%d', i);
    saveFigure(figure(i), fileName, 'fig');
end

%% Armazenamento dos resultados de simulação

save('results/HybridMPPT/test_case_out.mat', 'test_case_out', '-v7.3');
