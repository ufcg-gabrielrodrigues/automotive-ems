%% Partição da simulação

sim_split_flag = false;    	% Flag de particionamento

%% Parâmetros temporais

T_s = 1e-6; % Passo de cálculo utilizado pelo 'solver' [s]

%% Varredura de parâmetros

% Lista de parâmetros a serem varridos individualmente
i_f = 4.5;                      % Corrente de excitação [A]
T_list = [50 150 200]';         % Temperatura dos enrolamentos de estator [oC]
n_r_list = [2000 3500 7500]'; % Velocidade do alternador [rpm]
r_l_list = [0.5 1.0]';          % Resistência de carga [Ohm]

% Formação das casos de varredura
param_sweep = [];

for index_T = 1:length(T_list)
    
    n_r_sweep = [];
    [r_l_dim, ~] = size(r_l_list);
    
    for index_n_r = 1:length(n_r_list)
        n_r_tmp = n_r_list(index_n_r) * ones(r_l_dim, 1);
        n_r_tmp = [n_r_tmp, r_l_list];
        n_r_sweep = [n_r_sweep; n_r_tmp];
    end
    
    [n_r_sweep_dim, ~] = size(n_r_sweep);
    
    T_tmp = T_list(index_T) * ones(n_r_sweep_dim, 1);
    T_tmp = [T_tmp, n_r_sweep];
    param_sweep = [param_sweep; T_tmp];
end

%% Perfil temporal do ciclo de trabalho

% Período de simulação por valor de ciclo de trabalho
t_u = 1e-2; % [s]

% Formação da estrutura de dados contendo o perfil
u.Data = 0.0:0.01:1.0;
u.Time = [0.0 (t_u * 2):t_u:(t_u * length(u.Data))];

%% Inicializa modelo no Simulink

open_system('models/ThermalAnalysis.slx', 'loadonly');

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
set_param('ThermalAnalysis/Solver Configuration', 'UseLocalSolver', 'on');
set_param('ThermalAnalysis/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('ThermalAnalysis/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('ThermalAnalysis/Solver Configuration', 'DoFixedCost', 'on');
set_param('ThermalAnalysis/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
set_param('ThermalAnalysis', 'StopTime', num2str(u.Time(end) + t_u));

% Salva mundanças feitas no modelo
save_system('models/ThermalAnalysis.slx');

%% Configuração dos casos de teste como entrada do modelo no Simulink

[num_cases, ~] = size(param_sweep);

for test_case_index = 1:num_cases
    test_case = param_sweep(test_case_index, :);
    T = test_case(1);
    n_r = test_case(2);
    r_l = test_case(3);
    
    r_s = alternator.stator.r.function(T);
    
    simIn(test_case_index) = Simulink.SimulationInput('ThermalAnalysis');
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('ThermalAnalysis/i_f', 'Value', num2str(i_f));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('ThermalAnalysis/Alternator/Stator Resistance', 'constant', num2str(r_s));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('ThermalAnalysis/n_r', 'Value', num2str(n_r));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('ThermalAnalysis/r_l', 'R', num2str(r_l));
end

%% Execução da simulação em ambiente Simulink

simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
    'TransferBaseWorkspaceVariables', 'on');

%% Finaliza modelo no Simulink

close_system('models/ThermalAnalysis.slx');

%% Ajustes para iniciar tratamento e registro de resultados

% Ajuste de ciclo de trabalho pra inclusão de passo final
u.Time(end + 1) = u.Time(end) + t_u;
u.Data(end + 1) = u.Data(end);

%% Registro de resultados obtidos no caso de teste

% Laço de iterações por casos de teste
for test_case_index = 1:num_cases
    
    test_case = param_sweep(test_case_index, :);
    
    % Valor de variáveis correspondentes ao caso de teste
    T = test_case(1);
    n_r = test_case(2);
    r_l = test_case(3);
    
    % Alternador
    test_case_out(test_case_index).alternator.rotor.n = n_r;
    test_case_out(test_case_index).alternator.rotor.l.i = i_f;
    
    test_case_out(test_case_index).alternator.stator.r.T = T;
    test_case_out(test_case_index).alternator.stator.r.value = alternator.stator.r.function(T);
    test_case_out(test_case_index).alternator.stator.input.e.value = simOut(test_case_index).e_a_abc;
    test_case_out(test_case_index).alternator.stator.output.v = simOut(test_case_index).v_a_abc;
    test_case_out(test_case_index).alternator.stator.output.i = simOut(test_case_index).i_a_abc;
    
    % Retificador
    test_case_out(test_case_index).rectifier.control.u = timeseries();
    test_case_out(test_case_index).rectifier.control.u.Time = u.Time;
    test_case_out(test_case_index).rectifier.control.u.Data(:, 1) = u.Data;
    
    % Carga
    test_case_out(test_case_index).electrical_load.r = r_l;
    test_case_out(test_case_index).electrical_load.p.inst = simOut(test_case_index).p_l;
    test_case_out(test_case_index).electrical_load.p.avg = simOut(test_case_index).p_l;
    
    p_time = simOut(test_case_index).p_l.Time;
    
    for u_index = 2:length(u.Time)
        p_terms = find(p_time > (u.Time(u_index) - t_u) & p_time <= u.Time(u_index));
        p_piece = simOut(test_case_index).p_l.Data(p_terms);
        test_case_out(test_case_index).electrical_load.p.avg.Data(p_terms) = mean(p_piece);
    end
end

%% Tratamento de casos de teste

% Inicialização de variáveis pertinentes ao tratamento
t_interest = u.Time(1:(end - 1));
t_interest(1) = t_u;
t_interest = round((t_interest + t_u/2)/T_s);
case_length = length(u.Time) - 1;
test_case_matrix = zeros(num_cases * case_length, 5);

% Laço de iterações por casos de teste
for test_case_index = 1:num_cases
    
    case_interval = (case_length * (test_case_index - 1) + 1):(case_length * test_case_index);
    
    test_case_matrix(case_interval, 1) = param_sweep(test_case_index, 1);
    test_case_matrix(case_interval, 2) = test_case_out(test_case_index).alternator.rotor.n;
    test_case_matrix(case_interval, 3) = test_case_out(test_case_index).electrical_load.r;
    test_case_matrix(case_interval, 4) = test_case_out(test_case_index).rectifier.control.u.Data(1:(end - 1), 1);
    test_case_matrix(case_interval, 5) = test_case_out(test_case_index).electrical_load.p.avg.Data(t_interest, 1);
end

%% Identificação de pontos de máxima potência

mpp_u_3d = zeros(length(T_list), length(n_r_list), length(r_l_list));
mpp_p_3d = zeros(length(T_list), length(n_r_list), length(r_l_list));

mpp_matrix = [param_sweep zeros(num_cases, 2)];
case_index = 0;

T_index = 0;
n_r_index = 0;
r_l_index = 0;

for T = T_list'
    T_index = T_index + 1;
    n_r_index = 0;
    
    for n_r = n_r_list'
        n_r_index = n_r_index + 1;
        r_l_index = 0;
        
        for r_l = r_l_list'
            r_l_index = r_l_index + 1;
            case_index = case_index + 1;
            
            curve_indexes = find((test_case_matrix(:, 1) == T) ...
                & (test_case_matrix(:, 2) == n_r) ...
                & (test_case_matrix(:, 3) == r_l));
            
            u = test_case_matrix(curve_indexes, 4);
            p = test_case_matrix(curve_indexes, 5);
            
            [p_max, i_max] = max(p);
            u_max = u(i_max);
            
            mpp_matrix(case_index, 4) = u_max;
            mpp_matrix(case_index, 5) = p_max;
            
            mpp_u_3d(T_index, n_r_index, r_l_index) = u_max;
            mpp_p_3d(T_index, n_r_index, r_l_index) = p_max;
        end
        
    end
end

% % 
% [n_r, r_l] = meshgrid(n_r_list, r_l_list);
% mpp_u = reshape(mpp_u_3d, length(r_l_list), length(n_r_list));
% mpp_p = reshape(mpp_p_3d, length(r_l_list), length(n_r_list));
% 
% %
% figure_index = figure_index + 1;
% figure(figure_index)
% 
% surf(n_r, r_l, mpp_u);
% 
% % 
% figure_index = figure_index + 1;
% figure(figure_index)
% 
% surf(n_r, r_l, mpp_p);

%% Parâmetros auxiliares para figuras

figure_index = 0;
leg_entries_per_columns = 3;

%% Curvas de potência 

% Comparadas por temperatura de enrolamento do estator
for n_r = n_r_list'
    
    for r_l = r_l_list'
        
        colors = distinguishable_colors(length(T_list));
        color_index = 0;
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for T = T_list'
            color_index = color_index + 1;
            curve_indexes = find((test_case_matrix(:, 1) == T) ...
                & (test_case_matrix(:, 2) == n_r) ...
                & (test_case_matrix(:, 3) == r_l));
            
            plot(test_case_matrix(curve_indexes, 4), test_case_matrix(curve_indexes, 5), ...
                'Color', colors(color_index, :), 'DisplayName', ...
                ['$T = ' num2str(T) ' ^{o}C$']);
            hold on;
            
            legend('off');
            legend('show');
        end
        
        title(['Curvas de pot{\^{e}}ncia [$n_{r} = ' num2str(n_r) ...
            ' rpm$; $r_{l} = ' num2str(r_l) ' \Omega$]']);
        xlabel('$u$');
        ylabel('$P [W]$');
        leg = legend;
        leg.Location = 'SouthOutside';
        leg.NumColumns = ceil(length(T_list)/leg_entries_per_columns);
        grid on;
    end
end

%% Relação dos pontos de máxima potência com variáveis do sistema

% Relação com a temperatura de enrolamento do estator
for n_r = n_r_list'
    
    for r_l = r_l_list'
        
        curve_indexes = find((mpp_matrix(:, 2) == n_r) & (mpp_matrix(:, 3) == r_l));
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        subplot(2, 1, 1)
        plot(mpp_matrix(curve_indexes, 1), mpp_matrix(curve_indexes, 4), 'o-');
        xlabel('$T [^{o}C]$');
        ylabel('$u$');
        ylim([0 1]);
        grid on;
        
        subplot(2, 1, 2)
        plot(mpp_matrix(curve_indexes, 1), mpp_matrix(curve_indexes, 5), 'o-');
        xlabel('$T [^{o}C]$');
        ylabel('$P [W]$');
        grid on;
        
        suptitle(['Rela{\c{c}}{\~{a}}o entre o ponto de m{\''{a}}xima ' ...
            'pot{\^{e}}ncia e a temperatura de enrolamento do estator [$n_{r} = ' ...
            num2str(n_r) ' rpm$; $r_{l} = ' num2str(r_l) ' \Omega$]']);
    end
end

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = sprintf('results/ThermalAnalysis/Figura_%d', i);
    saveFigure(figure(i), fileName, 'fig');
end

%% Armazenamento dos resultados de simulação

save('results/ThermalAnalysis/test_case_out.mat', 'test_case_out', '-v7.3');
save('results/ThermalAnalysis/test_case_matrix.mat', 'test_case_matrix', '-v7.3');
save('results/ThermalAnalysis/mpp_matrix.mat', 'mpp_matrix', '-v7.3');
save('results/ThermalAnalysis/mpp_map.mat', 'T_list', 'n_r_list', 'r_l_list', ...
    'mpp_u_3d', 'mpp_p_3d', '-v7.3');
