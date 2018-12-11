%% Parâmetros temporais

T_s = 1e-6; % Passo de cálculo utilizado pelo 'solver' [s]

%% Retificador

% Filtro passivo
rectifier.filter.c = 10e-3;	% Capacitância de filtro [F]

%% Varredura de parâmetros

% Lista de parâmetros a serem varridos individualmente
i_f_list = 3.0;                             % Corrente de excitação [A]
n_alt_list = 2000:500:6000;                 % Velocidade do alternador [rpm]
r_l_list = [0.01 0.05:0.05:0.45 0.5:0.5:5];	% Resistência de carga [Ohm]

% Formação das casos de varredura
param_sweep = [];

for index_i_f = 1:length(i_f_list)
    
    n_alt_sweep = [];
    
    for index_n_alt = 1:length(n_alt_list)
        n_alt_tmp = n_alt_list(index_n_alt) * ones(length(r_l_list), 1);
        n_alt_tmp = [n_alt_tmp, r_l_list'];
        n_alt_sweep = [n_alt_sweep; n_alt_tmp];
    end
    
    i_f_tmp = i_f_list(index_i_f) * ones(length(n_alt_sweep), 1);
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

open_system('MPPTCurves.slx', 'loadonly');

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
set_param('MPPTCurves/Solver Configuration', 'UseLocalSolver', 'on');
set_param('MPPTCurves/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('MPPTCurves/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('MPPTCurves/Solver Configuration', 'DoFixedCost', 'on');
set_param('MPPTCurves/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
simulationParameters.StopTime = num2str(u.Time(end) + t_u); % [s]

%% Execução da simulação em ambiente Simulink

[num_cases, ~] = size(param_sweep);

for test_case_index = 1:num_cases
    test_case = param_sweep(test_case_index, :);
    i_f = test_case(1);
    n_alt = test_case(2);
    r_l = test_case(3);
    
    % 
    fprintf('Running test case %d/%d...\n', test_case_index, num_cases);
    
    %
    simout{test_case_index} = sim('MPPTCurves', simulationParameters);
end

%% Salva e finaliza modelo no Simulink

save_system('MPPTCurves.slx');
close_system('MPPTCurves.slx');

%% Registro de resultados obtidos no caso de teste

% 
u.Time(end + 1) = u.Time(end) + t_u;
u.Data(end + 1) = u.Data(end);

% 
for test_case_index = 1:num_cases
    % 
    test_case = param_sweep(test_case_index, :);
    i_f = test_case(1);
    n_alt = test_case(2);
    r_l = test_case(3);
    
    % Alternador
    test_case_out(test_case_index).alternator.rotor.n = n_alt;
    test_case_out(test_case_index).alternator.rotor.l.i = i_f;
    
    test_case_out(test_case_index).alternator.stator.input.e.value = simout{test_case_index}.e_a_abc;
    test_case_out(test_case_index).alternator.stator.output.v = simout{test_case_index}.v_a_abc;
    test_case_out(test_case_index).alternator.stator.output.i = simout{test_case_index}.i_a_abc;
    
    % Retificador
    test_case_out(test_case_index).rectifier.control.u = timeseries();
    test_case_out(test_case_index).rectifier.control.u.Time = u.Time;
    test_case_out(test_case_index).rectifier.control.u.Data(:, 1) = u.Data;
    
    % Carga
    test_case_out(test_case_index).load.r = r_l;
    test_case_out(test_case_index).load.p.inst = simout{test_case_index}.p_l;
    test_case_out(test_case_index).load.p.avg = simout{test_case_index}.p_l;
    
    p_time = simout{test_case_index}.p_l.Time;
    
    for u_index = 2:length(u.Time)
        p_terms = find(p_time > (u.Time(u_index) - t_u) & p_time <= u.Time(u_index));
        p_piece = simout{test_case_index}.p_l.Data(p_terms);
        test_case_out(test_case_index).load.p.avg.Data(p_terms) = mean(p_piece);
    end
end

%% Tratamento de casos de teste

t_interest = u.Time(1:(end - 1));
t_interest(1) = t_u;
t_interest = round((t_interest + t_u/2)/T_s);
batch_length = length(u.Time) - 1;
test_case_matrix = zeros(length(test_case_out) * batch_length, 5);

for test_case_index = 1:length(test_case_out)
    batch_interval = (batch_length * (test_case_index - 1) + 1):(batch_length * test_case_index);
    
    test_case_matrix(batch_interval, 1) = test_case_out(test_case_index).alternator.rotor.l.i;
    test_case_matrix(batch_interval, 2) = test_case_out(test_case_index).alternator.rotor.n;
    test_case_matrix(batch_interval, 3) = test_case_out(test_case_index).load.r;
    test_case_matrix(batch_interval, 4) = test_case_out(test_case_index).rectifier.control.u.Data(1:(end - 1), 1);
    test_case_matrix(batch_interval, 5) = test_case_out(test_case_index).load.p.avg.Data(t_interest, 1);
end

%% Parâmetros auxiliares para figuras

figure_index = 0;
leg_entries_per_columns = 3;

%% Curvas de potência (comparadas por velocidade do alternador)

for i_f = i_f_list
    
    for r_l = r_l_list
        
        colors = distinguishable_colors(length(n_alt_list));
        color_index = 0;
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for n_alt = n_alt_list
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
        ylabel('$P$ [W]');
        leg = legend;
        leg.Location = 'SouthOutside';
        leg.NumColumns = ceil(length(n_alt_list)/leg_entries_per_columns);
        grid on;
    end
end

%% Curvas de potência (comparadas por carga)

for i_f = i_f_list
    
    for n_alt = n_alt_list
        
        colors = distinguishable_colors(length(r_l_list));
        color_index = 0;
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for r_l = r_l_list
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
        ylabel('$P$ [W]');
        leg = legend;
        leg.Location = 'SouthOutside';
        leg.NumColumns = ceil(length(r_l_list)/leg_entries_per_columns);
        grid on;
    end
end

%% Identificação de pontos de máxima potência

mpp_u_3d = zeros(length(r_l_list), length(n_alt_list), length(i_f_list));
mpp_p_3d = zeros(length(r_l_list), length(n_alt_list), length(i_f_list));

mpp_matrix = [param_sweep zeros(num_cases, 2)];
case_index = 0;

i_f_index = 0;
n_alt_index = 0;
r_l_index = 0;

for i_f = i_f_list
    i_f_index = i_f_index + 1;
    n_alt_index = 0;
    
    for n_alt = n_alt_list
        n_alt_index = n_alt_index + 1;
        r_l_index = 0;
        
        for r_l = r_l_list
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
            
            mpp_u_3d(r_l_index, n_alt_index, i_f_index) = u_max;
            mpp_p_3d(r_l_index, n_alt_index, i_f_index) = p_max;
        end
        
    end
end

% 
[n_alt, r_l] = meshgrid(n_alt_list, r_l_list);
mpp_u = reshape(mpp_u_3d, length(r_l_list), length(n_alt_list));
mpp_p = reshape(mpp_p_3d, length(r_l_list), length(n_alt_list));

%
figure_index = figure_index + 1;
figure(figure_index)

surf(n_alt, r_l, mpp_u);

% 
figure_index = figure_index + 1;
figure(figure_index)

surf(n_alt, r_l, mpp_p);

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = sprintf('results/Figura_%d', i);
    saveFigure(figure(i), fileName, 'fig');
end

%% Armazenamento dos resultados de simulação

save('results/test_case_out.mat', 'test_case_out', '-v7.3');
save('results/test_case_matrix.mat', 'test_case_matrix', '-v7.3');
save('results/mpp_matrix.mat', 'mpp_matrix', '-v7.3');
