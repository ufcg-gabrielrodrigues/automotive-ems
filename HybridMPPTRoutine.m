%% Inicializa modelo no Simulink

open_system('models/HybridMPPT.slx', 'loadonly');

%% Parâmetros temporais

T_s = 1.0e-6;   % Passo de cálculo utilizado pelo 'solver' [s]
T_r = 1.0e-5;   % Passo de amostragem global de leitura de variáveis [s]
T_k = 1.0e-4;   % Passo de amostragem global de rotinas de controle [s]
t_f = 3.0e-1;   % Tempo total de simulação [s]

%% Esquema de controle

% Atualização de parâmetro: indutância mútua
m_f_default_local = 'm_f = 0;';

if (isfield(alternator.m_f, 'function'))
    m_f_str = regexprep(func2str(alternator.m_f.function), '@\(.+?\)', '');
else
    m_f_str = num2str(alternator.m_f.value);
end

if (alternator.stator.connection == delta)
    m_f_str = ['(' m_f_str ')./sqrt(3)'];
end

m_f_local = ['m_f = ' m_f_str ';'];
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'HybridMPPT/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, m_f_default_local, m_f_local);

%% Alternador

% Corrente de excitação máxima
i_f_max = 5;                    % [A]

% Temperatura da resistência do circuito de estator
alternator.stator.r.T = 150;    % [oC]

%% Retificador

% Parâmetros de controle
rectifier.control.pwm.f_s = 1/T_k;  % Frequência de chaveamento dos PWMs de controle do circuito retificador [Hz]

%% Bateria

battery.v_nom = 80.0;   % [V]

%% Varredura de parâmetros

% Lista de parâmetros a serem varridos individualmente
p_l_max = [4110; 5280; 6350; 7445];
n_r_list = (4500:1000:7500)';   % Velocidade do alternador [rpm]
k_p_list = [1]';                % 
k_i_list = [0 1]';              % 

% Formação das casos de varredura
param_sweep = combvec(n_r_list', k_p_list', k_i_list')';

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
    n_r = test_case(1);
    k_p = test_case(2);
    k_i = test_case(3);
    
    simIn(test_case_index) = Simulink.SimulationInput('HybridMPPT');
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('HybridMPPT/n_r', 'Value', num2str(n_r));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('HybridMPPT/k_p', 'Gain', num2str(k_p));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('HybridMPPT/k_i', 'Gain', num2str(k_i));
end

%% Execução da simulação em ambiente Simulink

% Execução da simulação paralelizada
simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
    'TransferBaseWorkspaceVariables', 'on');

%% Redefinição de parâmetros de alternador

% Atualização de parâmetro para valor padrão: indutância mútua
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'HybridMPPT/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, m_f_local, m_f_default_local);

%% Finaliza modelo no Simulink

save_system('models/HybridMPPT.slx');
close_system('models/HybridMPPT.slx');

%% Registro de resultados obtidos no caso de teste

% Laço de iterações por casos de teste
for test_case_index = 1:num_cases
    
    test_case = param_sweep(test_case_index, :);
    
    % Valor de variáveis correspondentes ao caso de teste
    n_r = test_case(1);
    k_p = test_case(2);
    k_i = test_case(3);
    
    % Esquema de controle
    test_case_out(test_case_index).k_p = k_p;
    test_case_out(test_case_index).k_i = k_i;
    
    % Alternador
    test_case_out(test_case_index).alternator.rotor.n = n_r;
    test_case_out(test_case_index).alternator.rotor.l.i = i_f_max;
    
    % Retificador
    test_case_out(test_case_index).rectifier.control.u = simOut(test_case_index).u_smr;
    
    % Carga
    test_case_out(test_case_index).electrical_load.v = simOut(test_case_index).v_l;
    test_case_out(test_case_index).electrical_load.i = simOut(test_case_index).i_l;
    test_case_out(test_case_index).electrical_load.p = simOut(test_case_index).p_l;
end

%% Verifica��o de an�lise de pot�ncia

if (length(p_l_max) == length(n_r_list))
    power_analysis_flag = true;
else
    power_analysis_flag = false;
end

%% Resultados comparativos

%
figure_index = 0;

%
cmp_dim = length(k_p_list)*length(k_i_list);

% Cores
if (power_analysis_flag)
    colors = distinguishable_colors(cmp_dim + 1);
else
    colors = distinguishable_colors(cmp_dim);
end

[color_dim, ~] = size(colors);

%
param_sweep = [param_sweep (1:length(param_sweep))'];
param_sweep = sortrows(param_sweep, 1);

for test_case_index = 1:num_cases/cmp_dim
    %
    test_case = param_sweep(1:cmp_dim, :);
    param_sweep(1:cmp_dim, :) = [];
    
    %
    n_r = test_case(1, 1);
    
    figure_index = figure_index + 1;
    figure(figure_index)
    
    subplot(2, 1, 1)
    
    for cmp_case = 1:cmp_dim
        k_p = test_case(cmp_case, 2);
        k_i = test_case(cmp_case, 3);
        
        if (k_p == 0 && k_i == 0)
            leg = '\textit{Load matching} desligado';
        elseif (k_p == 0 && k_i == 1)
            leg = 'ESC';
        elseif (k_p == 1 && k_i == 0)
            leg = 'Anal{\''{i}}tico';
        elseif (k_p == 1 && k_i == 1)
            leg = 'H{\''{i}}brido';
        end
            
        plot(test_case_out(test_case(cmp_case, end)).electrical_load.p.time, ...
            test_case_out(test_case(cmp_case, end)).electrical_load.p.data, ...
            'Color', colors(cmp_case, :), 'DisplayName', leg);
        
        legend('off');
        legend('show');
        
        hold on;
    end
    
    if (power_analysis_flag)
        plot([test_case_out(test_case(cmp_case, end)).electrical_load.p.time(1) ...
            test_case_out(test_case(cmp_case, end)).electrical_load.p.time(end)], ...
            p_l_max(test_case_index)*[1 1], '--', 'Color', colors(color_dim, :), ...
            'DisplayName', 'M{\''{a}}xima pot{\^{e}ncia}');
    end
    
    hold off;
    ylim(p_l_max(test_case_index)*[0.9 1.1]); 
    ylabel('$p_{l}\,[\textrm{W}]$');
    leg = legend;
    leg.Location = 'NorthEast';
    grid on;
    
    subplot(2, 1, 2)
    
    for cmp_case = 1:cmp_dim
        plot(test_case_out(test_case(cmp_case, end)).rectifier.control.u.time, ...
            test_case_out(test_case(cmp_case, end)).rectifier.control.u.data, ...
            'Color', colors(cmp_case, :));
        
        hold on;
    end
    
    hold off;
    ylim([0 Inf]);
    xlabel('$t\,[\textrm{s}]$');
    ylabel('$d_{\textrm{smr}}$');
    grid on;
    
    tit{test_case_index} = ['results/HybridMPPT/mppt-nr-' num2str(n_r)];
end

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = tit{i};
    saveFigure(figure(i), fileName, 'fig');
    saveFigure(figure(i), fileName, 'eps');
end

%% Armazenamento dos resultados de simulação

save('results/HybridMPPT/test_case_out.mat', 'test_case_out', '-v7.3');