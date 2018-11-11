%% Alternador

alternator.rotor.n = 2000;  % Valocidade do rotor [rpm]

%% Retificador

% Filtro passivo
rectifier.filter.c = 47e-3;	% Capacitância de filtro [F]

%% Carga

load.r = 13.5/35;          	% Resistência de carga [Ohm]

%% Inicializa modelo no Simulink

open_system('models/LundellAlternator.slx', 'loadonly');

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
T_s = 1e-6;
set_param('LundellAlternator/Solver Configuration', 'UseLocalSolver', 'on');
set_param('LundellAlternator/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('LundellAlternator/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('LundellAlternator/Solver Configuration', 'DoFixedCost', 'on');
set_param('LundellAlternator/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
simulationParameters.StopTime = '5e-1'; % [s]

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
load.v = simout.v_l;
load.i = simout.i_l;
load.p = simout.p_l;

%% Armazenamento dos resultados de simulação

save('results/alternator.mat', 'alternator', '-v7.3');
save('results/battery.mat', 'battery', '-v7.3');
save('results/load.mat', 'rectifier', '-v7.3');
