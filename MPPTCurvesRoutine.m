%% Partição da simulação

sim_split_flag = true;    	% Flag de particionamento
sim_batches = 60;           % Quantidade de partições

%% Armazenamento de resultados brutos

% Diretório
% raw_storage_path = 'results/MPPTCurves/';
raw_storage_path = 'D:/Users/gabriel.rodrigues/MATLAB/Projects/AutomotiveEMS/results/MPPTCurves/[19-01-15] Variant i_f - u_step 0.01/60-splitted/';

%% Parâmetros temporais

T_s = 1e-6; % Passo de cálculo utilizado pelo 'solver' [s]

%% Alternador

% Efeito t�rmico na resistência do circuito de estator
T = 150;                      	% [oC]
alternator.stator.r.value = alternator.stator.r.function(T);

%% Varredura de parâmetros

% Lista de parâmetros a serem varridos individualmente
i_f_list = (0.5:0.5:5.0)';                      % Corrente de excitação [A]
n_alt_list = (2000:500:7500)';                  % Velocidade do alternador [rpm]
r_l_list = [0.01 0.05:0.05:0.40 0.5:0.5:2.0]';	% Resistência de carga [Ohm]

% Formação das casos de varredura
param_sweep = [];

for index_i_f = 1:length(i_f_list)
    
    n_alt_sweep = [];
    [r_l_dim, ~] = size(r_l_list);
    
    for index_n_alt = 1:length(n_alt_list)
        n_alt_tmp = n_alt_list(index_n_alt) * ones(r_l_dim, 1);
        n_alt_tmp = [n_alt_tmp, r_l_list];
        n_alt_sweep = [n_alt_sweep; n_alt_tmp];
    end
    
    [n_alt_sweep_dim, ~] = size(n_alt_sweep);
    
    i_f_tmp = i_f_list(index_i_f) * ones(n_alt_sweep_dim, 1);
    i_f_tmp = [i_f_tmp, n_alt_sweep];
    param_sweep = [param_sweep; i_f_tmp];
end

%% Perfil temporal do ciclo de trabalho

% Período de simulação por valor de ciclo de trabalho
t_u = 2e-2; % [s]

% Formação da estrutura de dados contendo o perfil
u.Data = 0.0:0.01:1.0;
u.Time = [0.0 (t_u * 2):t_u:(t_u * length(u.Data))];

%% Inicializa modelo no Simulink

open_system('models/MPPTCurves.slx', 'loadonly');

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
set_param('MPPTCurves/Solver Configuration', 'UseLocalSolver', 'on');
set_param('MPPTCurves/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('MPPTCurves/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('MPPTCurves/Solver Configuration', 'DoFixedCost', 'on');
set_param('MPPTCurves/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
set_param('MPPTCurves', 'StopTime', num2str(u.Time(end) + t_u));

% Salva mundanças feitas no modelo
save_system('models/MPPTCurves.slx');

%% Configuração dos casos de teste como entrada do modelo no Simulink

[num_cases, ~] = size(param_sweep);

for test_case_index = 1:num_cases
    test_case = param_sweep(test_case_index, :);
    i_f = test_case(1);
    n_alt = test_case(2);
    r_l = test_case(3);
    
    simIn(test_case_index) = Simulink.SimulationInput('MPPTCurves');
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('MPPTCurves/i_f', 'Value', num2str(i_f));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('MPPTCurves/n_alt', 'Value', num2str(n_alt));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('MPPTCurves/r_l', 'R', num2str(r_l));
end

%% Execução da simulação em ambiente Simulink

% Verificação do particionamento da simulação
if (sim_split_flag)
    for batch_index = 1:sim_batches
        % Intervalo de simulação da iteração
        range = (num_cases/sim_batches*(batch_index - 1) + 1):(num_cases/sim_batches*batch_index);
        
        % Execução da simulação paralelizada
        simOut = parsim(simIn(range), 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
            'TransferBaseWorkspaceVariables', 'on');
        
        % Armazenamento dos resultados
        filename = [raw_storage_path 'simOut/simOut_' num2str(batch_index) '.mat'];
        save(filename, 'simOut', '-v7.3');
        
        % Limpeza das variáveis já armazenadas
        clear simOut;
    end
else
    % Redefinição do número de partições para uso subsequente
    sim_batches = 1;
    
    % Execução da simulação paralelizada
    simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
        'TransferBaseWorkspaceVariables', 'on');
end

%% Finaliza modelo no Simulink

close_system('models/MPPTCurves.slx');

%% Ajustes para iniciar tratamento e registro de resultados

% Ajuste de ciclo de trabalho pra inclusão de passo final
u.Time(end + 1) = u.Time(end) + t_u;
u.Data(end + 1) = u.Data(end);

%% Registro de resultados obtidos no caso de teste

% Laço de iterações por partição
for batch_index = 1:sim_batches
    
    % Tratamento para simulação particionada
    if (sim_split_flag)
        filename = [raw_storage_path 'simOut/simOut_' num2str(batch_index) '.mat'];
        load(filename);
    end
    
    % Laço de iterações por casos de teste
    for test_case_index = 1:num_cases/sim_batches
        
        % Tratamento para simulação particionada
        if (sim_split_flag)
            composed_test_case_index = (batch_index - 1) * num_cases/sim_batches + test_case_index;
            test_case = param_sweep(composed_test_case_index, :);
        else
            test_case = param_sweep(test_case_index, :);
        end
        
        % Valor de variáveis correspondentes ao caso de teste
        i_f = test_case(1);
        n_alt = test_case(2);
        r_l = test_case(3);
        
        % Alternador
        test_case_out(test_case_index).alternator.rotor.n = n_alt;
        test_case_out(test_case_index).alternator.rotor.l.i = i_f;
        
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
    
    % Tratamento para simulação particionada
    if (sim_split_flag)
        clear simOut;
        
        filename = [raw_storage_path 'test_case_out/test_case_out_' num2str(batch_index) '.mat'];
        save(filename, 'test_case_out', '-v7.3');
        clear test_case_out;
    end
end

%% Tratamento de casos de teste

% Inicialização de variáveis pertinentes ao tratamento
t_interest = u.Time(1:(end - 1));
t_interest(1) = t_u;
t_interest = round((t_interest + t_u/2)/T_s);
case_length = length(u.Time) - 1;
test_case_matrix = zeros(num_cases * case_length, 5);

% Laço de iterações por partição
for batch_index = 1:sim_batches
    
    % Tratamento para simulação particionada
    if (sim_split_flag)
        filename = [raw_storage_path 'test_case_out/test_case_out_' num2str(batch_index) '.mat'];
        load(filename);
    end
    
    % Laço de iterações por casos de teste
    for test_case_index = 1:num_cases/sim_batches
        
        % Tratamento para simulação particionada
        if (sim_split_flag)
            composed_test_case_index = (batch_index - 1) * num_cases/sim_batches + test_case_index;
            case_interval = (case_length * (composed_test_case_index - 1) + 1):(case_length * composed_test_case_index);
        else
            case_interval = (case_length * (test_case_index - 1) + 1):(case_length * test_case_index);
        end
        
        test_case_matrix(case_interval, 1) = test_case_out(test_case_index).alternator.rotor.l.i;
        test_case_matrix(case_interval, 2) = test_case_out(test_case_index).alternator.rotor.n;
        test_case_matrix(case_interval, 3) = test_case_out(test_case_index).electrical_load.r;
        test_case_matrix(case_interval, 4) = test_case_out(test_case_index).rectifier.control.u.Data(1:(end - 1), 1);
        test_case_matrix(case_interval, 5) = test_case_out(test_case_index).electrical_load.p.avg.Data(t_interest, 1);
    end
    
    % Tratamento para simulação particionada
    if (sim_split_flag)
        clear test_case_out;
    end
end

%% Identificação de pontos de máxima potência

mpp_u_3d = zeros(length(i_f_list), length(n_alt_list), length(r_l_list));
mpp_p_3d = zeros(length(i_f_list), length(n_alt_list), length(r_l_list));

mpp_matrix = [param_sweep zeros(num_cases, 2)];
case_index = 0;

i_f_index = 0;
n_alt_index = 0;
r_l_index = 0;

for i_f = i_f_list'
    i_f_index = i_f_index + 1;
    n_alt_index = 0;
    
    for n_alt = n_alt_list'
        n_alt_index = n_alt_index + 1;
        r_l_index = 0;
        
        for r_l = r_l_list'
            r_l_index = r_l_index + 1;
            case_index = case_index + 1;
            
            curve_indexes = find((test_case_matrix(:, 1) == i_f) ...
                & (test_case_matrix(:, 2) == n_alt) ...
                & (test_case_matrix(:, 3) == r_l));
            
            u = test_case_matrix(curve_indexes, 4);
            p = test_case_matrix(curve_indexes, 5);
            
            [p_max, i_max] = max(p);
            u_max = u(i_max);
            
            mpp_matrix(case_index, 4) = u_max;
            mpp_matrix(case_index, 5) = p_max;
            
            mpp_u_3d(i_f_index, n_alt_index, r_l_index) = u_max;
            mpp_p_3d(i_f_index, n_alt_index, r_l_index) = p_max;
        end
        
    end
end

% % 
% [n_alt, r_l] = meshgrid(n_alt_list, r_l_list);
% mpp_u = reshape(mpp_u_3d, length(r_l_list), length(n_alt_list));
% mpp_p = reshape(mpp_p_3d, length(r_l_list), length(n_alt_list));
% 
% %
% figure_index = figure_index + 1;
% figure(figure_index)
% 
% surf(n_alt, r_l, mpp_u);
% 
% % 
% figure_index = figure_index + 1;
% figure(figure_index)
% 
% surf(n_alt, r_l, mpp_p);

%% Parâmetros auxiliares para figuras

figure_index = 0;
leg_entries_per_columns = 3;

%% Curvas de potência 

% Comparadas por corrente de excitação
for n_alt = n_alt_list'
    
    for r_l = r_l_list'
        
        colors = distinguishable_colors(length(i_f_list));
        color_index = 0;
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for i_f = i_f_list'
            color_index = color_index + 1;
            curve_indexes = find((test_case_matrix(:, 1) == i_f) ...
                & (test_case_matrix(:, 2) == n_alt) ...
                & (test_case_matrix(:, 3) == r_l));
            
            plot(test_case_matrix(curve_indexes, 4), test_case_matrix(curve_indexes, 5), ...
                'Color', colors(color_index, :), 'DisplayName', ...
                ['$i_{f} = ' num2str(i_f) ' A$']);
            hold on;
            
            legend('off');
            legend('show');
        end
        
        title(['Curvas de pot{\^{e}}ncia [$n_{alt} = ' num2str(n_alt) ...
            ' rpm$; $r_{l} = ' num2str(r_l) ' \Omega$]']);
        xlabel('$u$');
        ylabel('$P [W]$');
        leg = legend;
        leg.Location = 'SouthOutside';
        leg.NumColumns = ceil(length(i_f_list)/leg_entries_per_columns);
        grid on;
    end
end

% Comparadas por velocidade do alternador
for i_f = i_f_list'
    
    for r_l = r_l_list'
        
        colors = distinguishable_colors(length(n_alt_list));
        color_index = 0;
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for n_alt = n_alt_list'
            color_index = color_index + 1;
            curve_indexes = find((test_case_matrix(:, 1) == i_f) ...
                & (test_case_matrix(:, 2) == n_alt) ...
                & (test_case_matrix(:, 3) == r_l));
            
            plot(test_case_matrix(curve_indexes, 4), test_case_matrix(curve_indexes, 5), ...
                'Color', colors(color_index, :), 'DisplayName', ...
                ['$n_{alt} = ' num2str(n_alt) ' rpm$']);
            hold on;
            
            legend('off');
            legend('show');
        end
        
        title(['Curvas de pot{\^{e}}ncia [$i_{f} = ' num2str(i_f) ...
            ' A$; $r_{l} = ' num2str(r_l) ' \Omega$]']);
        xlabel('$u$');
        ylabel('$P [W]$');
        leg = legend;
        leg.Location = 'SouthOutside';
        leg.NumColumns = ceil(length(n_alt_list)/leg_entries_per_columns);
        grid on;
    end
end

% Comparadas por carga
for i_f = i_f_list'
    
    for n_alt = n_alt_list'
        
        colors = distinguishable_colors(length(r_l_list));
        color_index = 0;
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for r_l = r_l_list'
            color_index = color_index + 1;
            curve_indexes = find((test_case_matrix(:, 1) == i_f) ...
                & (test_case_matrix(:, 2) == n_alt) ...
                & (test_case_matrix(:, 3) == r_l));
            
            plot(test_case_matrix(curve_indexes, 4), test_case_matrix(curve_indexes, 5), ...
                'Color', colors(color_index, :), 'DisplayName', ...
                ['$r_{l} = ' num2str(r_l, '%1.2f') ' \Omega$']);
            hold on;
            
            legend('off');
            legend('show');
        end
        
        title(['Curvas de pot{\^{e}}ncia [$i_{f} = ' num2str(i_f) ...
            ' A$; $n_{alt} = ' num2str(n_alt) ' rpm$]']);
        xlabel('$u$');
        ylabel('$P [W]$');
        leg = legend;
        leg.Location = 'SouthOutside';
        leg.NumColumns = ceil(length(r_l_list)/leg_entries_per_columns);
        grid on;
    end
end

%% Relação dos pontos de máxima potência com variáveis do sistema

% Relação com a corrente de excitação
for n_alt = n_alt_list'
    
    for r_l = r_l_list'
        
        curve_indexes = find((mpp_matrix(:, 2) == n_alt) & (mpp_matrix(:, 3) == r_l));
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        subplot(2, 1, 1)
        plot(mpp_matrix(curve_indexes, 1), mpp_matrix(curve_indexes, 4), 'o-');
        xlabel('$i_{f} [A]$');
        ylabel('$u$');
        ylim([0 1]);
        grid on;
        
        subplot(2, 1, 2)
        plot(mpp_matrix(curve_indexes, 1), mpp_matrix(curve_indexes, 5), 'o-');
        xlabel('$i_{f} [A]$');
        ylabel('$P [W]$');
        grid on;
        
        suptitle(['Rela{\c{c}}{\~{a}}o entre o ponto de m{\''{a}}xima ' ...
            'pot{\^{e}}ncia e a corrente de excita{\c{c}}{\~{a}}o [$n_{alt} = ' ...
            num2str(n_alt) ' rpm$; $r_{l} = ' num2str(r_l) ' \Omega$]']);
    end
end

% Relação com a velocidade do alternador
for i_f = i_f_list'
    
    for r_l = r_l_list'
        
        curve_indexes = find((mpp_matrix(:, 1) == i_f) & (mpp_matrix(:, 3) == r_l));
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        subplot(2, 1, 1)
        plot(mpp_matrix(curve_indexes, 2), mpp_matrix(curve_indexes, 4), 'o-');
        xlabel('$n_{alt} [rpm]$');
        ylabel('$u$');
        ylim([0 1]);
        grid on;
        
        subplot(2, 1, 2)
        plot(mpp_matrix(curve_indexes, 2), mpp_matrix(curve_indexes, 5), 'o-');
        xlabel('$n_{alt} [rpm]$');
        ylabel('$P [W]$');
        grid on;
        
        suptitle(['Rela{\c{c}}{\~{a}}o entre o ponto de m{\''{a}}xima ' ...
            'pot{\^{e}}ncia e a velocidade do alternador [$i_{f} = ' ...
            num2str(i_f) ' A$; $r_{l} = ' num2str(r_l) ' \Omega$]']);
    end
end

% Relação com a impedância de carga
for i_f = i_f_list'
    
    for n_alt = n_alt_list'
        
        curve_indexes = find((mpp_matrix(:, 1) == i_f) & (mpp_matrix(:, 2) == n_alt));
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        subplot(2, 1, 1)
        plot(mpp_matrix(curve_indexes, 3), mpp_matrix(curve_indexes, 4), 'o-');
        xlabel('$r_{l} [\Omega]$');
        ylabel('$u$');
        ylim([0 1]);
        grid on;
        
        subplot(2, 1, 2)
        plot(mpp_matrix(curve_indexes, 3), mpp_matrix(curve_indexes, 5), 'o-');
        xlabel('$r_{l} [\Omega]$');
        ylabel('$P [W]$');
        grid on;
        
        suptitle(['Rela{\c{c}}{\~{a}}o entre o ponto de m{\''{a}}xima ' ...
            'pot{\^{e}}ncia e a imped{\^{a}}ncia de carga [$i_{f} = ' ...
            num2str(i_f) ' A$; $n_{alt} = ' num2str(n_alt) ' rpm$]']);
    end
end

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = sprintf('results/MPPTCurves/Figura_%d', i);
    saveFigure(figure(i), fileName, 'fig');
end

%% Armazenamento dos resultados de simulação

% Tratamento para simulação particionada
if (~sim_split_flag)
    filename = [raw_storage_path 'test_case_out.mat'];
    save(filename, 'test_case_out', '-v7.3');
end

save('results/MPPTCurves/test_case_matrix.mat', 'test_case_matrix', '-v7.3');
save('results/MPPTCurves/mpp_matrix.mat', 'mpp_matrix', '-v7.3');
save('results/MPPTCurves/mpp_map.mat', 'i_f_list', 'n_alt_list', 'r_l_list', ...
    'mpp_u_3d', 'mpp_p_3d', '-v7.3');
