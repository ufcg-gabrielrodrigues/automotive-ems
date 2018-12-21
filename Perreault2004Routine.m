%% Inicializa modelo no Simulink

open_system('models/Perreault2004.slx', 'loadonly');

%% Parâmetros temporais

T_s = 1e-6; % Passo de cálculo utilizado pelo 'solver' local para sistemas físicos [s]
t_f = 0.1;  % Tempo total de simulação [s]

%% Motor a combustão interna

% Razão entre a polia do motor a combustão e a polia do alternador
iceToAltRotRatio = 2.5;

% Velocidade
n_ice = 1800/iceToAltRotRatio;

%% Alternador

% Corrente de excita��o m�xima
i_f_max = 3.6;  % [A]

% Atualização de parâmetro: fator de acoplamento
k_v_str = regexprep(func2str(alternator.k_v), '@\(.+?\)', '');
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'Perreault2004/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, 'k_v = 0', ['k_v = ' k_v_str]);

%% Retificador

% Filtro passivo
rectifier.filter.c = 10e-3;	% Capacitância de filtro [F]

%% Carga el�trica

electrical_load.r = 0.15;   % [Ohm]

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
set_param('Perreault2004/Solver Configuration', 'UseLocalSolver', 'on');
set_param('Perreault2004/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('Perreault2004/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('Perreault2004/Solver Configuration', 'DoFixedCost', 'on');
set_param('Perreault2004/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
simulationParameters.StopTime = num2str(t_f);   % [s]

%% Execução da simulação em ambiente Simulink

simout = sim('Perreault2004', simulationParameters);

%% Redefinição de parâmetros de alternador

% Atualização de parâmetro para valor padrão: fator de acoplamento
blockHandle.Script = strrep(blockHandle.Script, ['k_v = ' k_v_str], 'k_v = 0');

%% Salva e finaliza modelo no Simulink

save_system('models/Perreault2004.slx');
close_system('models/Perreault2004.slx');

%% Registro de resultados obtidos na simulação

% Motor a combustão interna
ice.n = n_ice;

% Alternador
alternator.rotor.n = simout.n_r;
alternator.rotor.l.i = i_f_max;

alternator.stator.input.e.value = simout.e_a_abc;
alternator.stator.output.v = simout.v_a_abc;
alternator.stator.output.i = simout.i_a_abc;

% Retificador
rectifier.control.u = simout.u_smr;

% Carga
electrical_load.v = simout.v_l;
electrical_load.i = simout.i_l;
electrical_load.z = simout.z_l;
electrical_load.p = simout.p_l;

%% Armazenamento dos resultados de simulação

save('results/Perreault2004/ice.mat', 'ice', '-v7.3');
save('results/Perreault2004/alternator.mat', 'alternator', '-v7.3');
save('results/Perreault2004/rectifier.mat', 'rectifier', '-v7.3');
save('results/Perreault2004/electrical_load.mat', 'electrical_load', '-v7.3');
