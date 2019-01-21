%% Inicializa modelo no Simulink

open_system('models/Perreault2004.slx', 'loadonly');

%% Parâmetros temporais

T_s = 1.0e-7;   % Passo de cálculo utilizado pelo 'solver' local para sistemas físicos [s]
T_k = 1.0e-5;   % Passo de amostragem global de rotinas de controle [s]
t_f = 2.0e-2;   % Tempo total de simulação [s]

%% Alternador

% Velocidade
alternator.rotor.n = 5000;

% Corrente de excita��o m�xima
alternator.rotor.l.i = 3.6; % [A]

% Atualização de parâmetro: fator de acoplamento
if (isfield(alternator.k_e, 'function'))
    k_e_str = regexprep(func2str(alternator.k_e.function), '@\(.+?\)', '');
else
    k_e_str = num2str(alternator.k_e.value);
end
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'Perreault2004/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, 'k_e = 0;', ['k_e = ' k_e_str ';']);

%% Carga el�trica

electrical_load.battery.v_nom = 50.0;   % [V]
electrical_load.r = 1.00;              	% [Ohm]

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
blockHandle.Script = strrep(blockHandle.Script, ['k_e = ' k_e_str ';'], 'k_e = 0;');

%% Salva e finaliza modelo no Simulink

save_system('models/Perreault2004.slx');
close_system('models/Perreault2004.slx');

%% Registro de resultados obtidos na simulação

% Alternador
alternator.stator.input.e = simout.e_a_abc;
alternator.stator.output.v = simout.v_a_abc;
alternator.stator.output.i = simout.i_a_abc;

% Retificador
rectifier.control.u = simout.u_smr;

% Carga
electrical_load.v = simout.v_l;
electrical_load.i = simout.i_l;
electrical_load.p = simout.p_l;

%% Armazenamento dos resultados de simulação

save('results/Perreault2004/alternator.mat', 'alternator', '-v7.3');
save('results/Perreault2004/rectifier.mat', 'rectifier', '-v7.3');
save('results/Perreault2004/electrical_load.mat', 'electrical_load', '-v7.3');
