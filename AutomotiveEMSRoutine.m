%% 

n_ice_option = 3;
iceToAltRotRatio = 2.5;

%% Inicializa modelo no Simulink

open_system('models/AutomotiveEMS.slx', 'loadonly');

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
T_s = 1e-6;
set_param('AutomotiveEMS/Solver Configuration', 'UseLocalSolver', 'on');
set_param('AutomotiveEMS/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('AutomotiveEMS/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('AutomotiveEMS/Solver Configuration', 'DoFixedCost', 'on');
set_param('AutomotiveEMS/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
simulationParameters.StopTime = '3e+0'; % [s]

%% Execução da simulação em ambiente Simulink

simout = sim('AutomotiveEMS', simulationParameters);

%% Salva e finaliza modelo no Simulink

save_system('models/AutomotiveEMS.slx');
close_system('models/AutomotiveEMS.slx');

%% Registro de resultados obtidos na simulação

% Motor a combustão interna
ice.n = simout.n_ice;

% Alternador
alternator.rotor.control.u = simout.u_i_f;
alternator.rotor.l.i = simout.i_f;

alternator.stator.input.e.value = simout.e_a_abc;
alternator.stator.output.v = simout.v_a_abc;
alternator.stator.output.i = simout.i_a_abc;

% Retificador
rectifier.control.u = simout.u_v_r;
rectifier.output.v = timeseries(simout.rectifier_output.Data(:, 1), simout.rectifier_output.Time);
rectifier.output.i = timeseries(simout.rectifier_output.Data(:, 2), simout.rectifier_output.Time);

% Bateria
battery.v = simout.v_b;

%% Armazenamento dos resultados de simulação

save('results/ice.mat', 'ice', '-v7.3');
save('results/alternator.mat', 'alternator', '-v7.3');
save('results/rectifier.mat', 'rectifier', '-v7.3');
save('results/battery.mat', 'battery', '-v7.3');
