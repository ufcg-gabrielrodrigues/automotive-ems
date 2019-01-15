%% Parâmetros temporais

T_s = 1e-6; % Passo de cálculo utilizado pelo 'solver' [s]
t_f = 5e-1; % Tempo total de simulação [s]

%% Alternador

% Corrente de excita��o m�xima
i_f_max = 5.0e-0;               % [A]

% Efeito t�rmico na resist�ncia do circuito de estator
T = 32;                         % [oC]
alternator.stator.r.value = alternator.stator.r.function(T);

%% Inicializa modelo no Simulink

open_system('models/LundellAlternator.slx', 'loadonly');

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
set_param('LundellAlternator/Solver Configuration', 'UseLocalSolver', 'on');
set_param('LundellAlternator/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('LundellAlternator/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('LundellAlternator/Solver Configuration', 'DoFixedCost', 'on');
set_param('LundellAlternator/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
set_param('LundellAlternator', 'StopTime', num2str(t_f));

% Salva mundanças feitas no modelo
save_system('models/LundellAlternator.slx');

%% Configuração dos casos de teste como entrada do modelo no Simulink

% Varredura de parâmetros [n_r i_dc v_dc]
param_sweep = [1967 36.058 13.49;
               3994 36.240 13.54;
               3983 55.650 13.58;
               5992 36.220 13.61;
               5985 56.190 13.46;
               5967 75.610 13.58];

% 
[num_cases, ~] = size(param_sweep);

for test_case_index = 1:num_cases
    % 
    test_case = param_sweep(test_case_index, :);
    n_r = test_case(1);
    i_dc = test_case(2);
    v_dc = test_case(3);
    
    % 
    r_dc = v_dc/i_dc;
    
    simIn(test_case_index) = Simulink.SimulationInput('LundellAlternator');
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('LundellAlternator/n_r', 'Value', num2str(n_r));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('LundellAlternator/r_l', 'R', num2str(r_dc));
end

%% Execução da simulação em ambiente Simulink

% Execução da simulação paralelizada
simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
    'TransferBaseWorkspaceVariables', 'on');

%% Salva e finaliza modelo no Simulink

save_system('models/LundellAlternator.slx');
close_system('models/LundellAlternator.slx');

%% Registro de resultados obtidos na simulação

% Laço de iterações por casos de teste
for test_case_index = 1:num_cases
    % Alternador
    test_case_out(test_case_index).alternator.rotor.n = simOut(test_case_index).n_r;
    test_case_out(test_case_index).alternator.rotor.l.i = simOut(test_case_index).i_f;

    test_case_out(test_case_index).alternator.stator.input.e = simOut(test_case_index).e_a_abc;
    test_case_out(test_case_index).alternator.stator.output.v = simOut(test_case_index).v_a_abc;
    test_case_out(test_case_index).alternator.stator.output.i = simOut(test_case_index).i_a_abc;
    test_case_out(test_case_index).alternator.stator.output.v_ll = simOut(test_case_index).v_ll_rms;
    test_case_out(test_case_index).alternator.stator.output.i_l = simOut(test_case_index).i_l_rms;

    % Bateria
    test_case_out(test_case_index).battery.v = simOut(test_case_index).v_b;

    % Carga
    test_case_out(test_case_index).electrical_load.v = simOut(test_case_index).v_l;
    test_case_out(test_case_index).electrical_load.i = simOut(test_case_index).i_l;
    test_case_out(test_case_index).electrical_load.p = simOut(test_case_index).p_l;
end

%% Armazenamento dos resultados de simulação

save('results/LundellAlternator/test_case_out.mat', 'test_case_out', '-v7.3');
