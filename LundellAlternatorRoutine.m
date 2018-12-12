%% Parâmetros temporais

T_s = 1e-6; % Passo de cálculo utilizado pelo 'solver' [s]
t_f = 5e-1; % Tempo total de simulação [s]

%% Alternador

alternator.rotor.n = 2000;  % Valocidade do rotor [rpm]

%% Retificador

% Filtro passivo
rectifier.filter.c = 47e-3;	% Capacitância de filtro [F]

%% Carga el�trica

electrical_load.r = 13.5/35;          	% Resistência de carga [Ohm]

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
simulationParameters.StopTime = num2str(t_f);   % [s]

%% Execução da simulação em ambiente Simulink

simout = sim('LundellAlternator', simulationParameters);

%% Salva e finaliza modelo no Simulink

save_system('models/LundellAlternator.slx');
close_system('models/LundellAlternator.slx');

%% Registro de resultados obtidos na simulação

% Alternador
alternator.rotor.n = simout.n_r;
alternator.rotor.l.i = simout.i_f;

alternator.stator.input.e.value = simout.e_a_abc;
alternator.stator.output.v = simout.v_a_abc;
alternator.stator.output.i = simout.i_a_abc;

% Bateria
battery.v = simout.v_b;

% Carga
electrical_load.v = simout.v_l;
electrical_load.i = simout.i_l;
electrical_load.p = simout.p_l;

%% Armazenamento dos resultados de simulação

save('results/alternator.mat', 'alternator', '-v7.3');
save('results/rectifier.mat', 'rectifier', '-v7.3');
save('results/electrical_load.mat', 'electrical_load', '-v7.3');
save('results/battery.mat', 'battery', '-v7.3');
