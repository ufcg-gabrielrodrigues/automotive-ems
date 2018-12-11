%% Parâmetros temporais

T_s = 1e-6; % Passo de cálculo utilizado pelo 'solver' [s]
t_f = 5.0;  % Tempo total de simula��o [s]

%% Motor a combustão interna

% Razão entre a polia do motor a combustão e a polia do alternador
iceToAltRotRatio = 2.5;

% Pontos de interesse do perfil de velocidade
n_ice_i = 6e+3/iceToAltRotRatio;
n_ice_f = 2e+3/iceToAltRotRatio;
t_brake_i = t_f*0.20;
t_brake_f = t_f*0.80;

%% Retificador

% Filtro passivo
rectifier.filter.c = 10e-3;	% Capacitância de filtro [F]

%% Carga

load.r = 0.5;

%% Esquemas de controle

conventionalControlScheme = Simulink.Variant('control_scheme == 0');
neuralControlScheme = Simulink.Variant('control_scheme == 1');
fittedControlScheme = Simulink.Variant('control_scheme == 2');

%% Inicializa modelo no Simulink

open_system('models/ControlComparison.slx', 'loadonly');

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
set_param('ControlComparison/Solver Configuration', 'UseLocalSolver', 'on');
set_param('ControlComparison/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('ControlComparison/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('ControlComparison/Solver Configuration', 'DoFixedCost', 'on');
set_param('ControlComparison/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
simulationParameters.StopTime = num2str(t_f);   % [s]

%% Execução da simulação em ambiente Simulink

parfor control_scheme_index = 1:3
    control_scheme = control_scheme_index;
    simout{control_scheme_index} = sim('ControlComparison', simulationParameters);
end

%% Salva e finaliza modelo no Simulink

save_system('models/ControlComparison.slx');
close_system('models/ControlComparison.slx');

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
rectifier.control.u = simout.u_smr;
rectifier.output.v = simout.v_r_o;
rectifier.output.i = simout.i_r_o;

% Bateria
battery.v = simout.v_b;

%% Armazenamento dos resultados de simulação

save('results/ice.mat', 'ice', '-v7.3');
save('results/alternator.mat', 'alternator', '-v7.3');
save('results/rectifier.mat', 'rectifier', '-v7.3');
save('results/battery.mat', 'battery', '-v7.3');
