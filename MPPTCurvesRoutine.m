 %% Retificador

% Filtro passivo
rectifier.filter.c = 47e-3;	% Capacitância de filtro [F]

%% Varredura de parâmetros

% Lista de parâmetros a serem varridos individualmente
i_f_list = 3.5;                 % Corrente de excitação [A]
n_alt_list = 2000:500:6000;     % Velocidade do alternador [rpm]
r_l_list = [0.01 0.5:0.5:5];	% Resistência de carga [Ohm]

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
T_s = 1e-6;
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

%% Armazenamento dos resultados de simulação

save('results/test_case_out.mat', 'test_case_out', '-v7.3');
save('results/test_case_matrix.mat', 'test_case_matrix', '-v7.3');
