%% Inicializa modelo no Simulink

open_system('models/AutomotiveEMS.slx', 'loadonly');

%% Parâmetros temporais

T_s = 1.0e-6;   % Passo de cálculo utilizado pelo 'solver' [s]
T_r = 1.0e-5;   % Passo de amostragem global de leitura de variáveis [s]
T_k = 1.0e-4;   % Passo de amostragem global de rotinas de controle [s]
t_f = 3.1e-0;   % Tempo total de simulação [s]

%% Motor a combustão interna

% Razão entre a polia do motor a combustão e a polia do alternador
iceToAltRotRatio = 2.5;

% Pontos de interesse do perfil de velocidade
n_ice_i = 7.5e+3/iceToAltRotRatio;
n_ice_f = 2.0e+3/iceToAltRotRatio;
t_brake = 3.0e-0;                   % [s]
t_brake_i = (t_f - t_brake);        % [s]
t_brake_f = t_brake_i + t_brake;    % [s]

%% Esquema de controle

% Atualização de parâmetro: fator de acoplamento
k_e_default_local = 'k_e = 0;';

if (isfield(alternator.k_e, 'function'))
    k_e_str = regexprep(func2str(alternator.k_e.function), '@\(.+?\)', '');
else
    k_e_str = num2str(alternator.k_e.value);
end

if (alternator.stator.connection == delta)
    k_e_str = ['(' k_e_str ')./sqrt(3)'];
end

k_e_local = ['k_e = ' k_e_str ';'];
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'AutomotiveEMS/Hybrid Load Matching Controller Scheme/Hybrid Load Matching Controller/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, k_e_default_local, k_e_local);

%% Alternador

% Corrente de excitação máxima
i_f_max = 4.75;                 % [A]

% Temperatura da resistência do circuito de estator
alternator.stator.r.T = 150;    % [oC]

%% Retificador

% Parâmetros de controle
rectifier.control.pwm.f_s = 1/T_k;  % Frequência de chaveamento dos PWMs de controle do circuito retificador [Hz]

%% Bateria

battery.v_nom = 12.0;   % [V]
battery.r = 50e-3;      % [Ohm]

%% Carga elétrica

electrical_load.r = 1.0e-0;	% [Ohm]

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
set_param('AutomotiveEMS/Solver Configuration', 'UseLocalSolver', 'on');
set_param('AutomotiveEMS/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('AutomotiveEMS/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('AutomotiveEMS/Solver Configuration', 'DoFixedCost', 'on');
set_param('AutomotiveEMS/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
simulationParameters.StopTime = num2str(t_f);   % [s]

% Salva mundanças feitas no modelo
save_system('models/AutomotiveEMS.slx');

%% Execução da simulação em ambiente Simulink

simOut = sim('AutomotiveEMS', simulationParameters);

%% Redefinição de parâmetros de alternador

% Atualização de parâmetro para valor padrão: fator de acoplamento
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'AutomotiveEMS/Hybrid Load Matching Controller Scheme/Hybrid Load Matching Controller/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, k_e_local, k_e_default_local);

%% Salva e finaliza modelo no Simulink

save_system('models/AutomotiveEMS.slx');
close_system('models/AutomotiveEMS.slx');

%% Registro de resultados obtidos na simulação

% Motor a combustão interna
ice.n = simOut.n_ice;

% Alternador
alternator.rotor.control.u = simOut.u_i_f;
alternator.rotor.l.i = simOut.i_f;

alternator.stator.input.e = simOut.e_a_abc;
alternator.stator.output.v = simOut.v_a_abc;
alternator.stator.output.i = simOut.i_a_abc;

% Retificador
rectifier.control.u = simOut.u_v_r;

% Carga elétrica
electrical_load.v = simOut.v_l;
electrical_load.i = simOut.i_l;
electrical_load.p = simOut.p_l;

% Barramento de circuito secundário
bus.v = simOut.v_bus;
bus.i = simOut.i_bus;
bus.p = simOut.p_bus;

% Banco de supercapacitores
uc_bank.v = simOut.v_uc;
uc_bank.i = simOut.i_uc;
uc_bank.p = simOut.p_uc;

%% Armazenamento dos resultados de simulação

save('results/AutomotiveEMS/ice.mat', 'ice', '-v7.3');
save('results/AutomotiveEMS/alternator.mat', 'alternator', '-v7.3');
save('results/AutomotiveEMS/rectifier.mat', 'rectifier', '-v7.3');
save('results/AutomotiveEMS/battery.mat', 'battery', '-v7.3');
save('results/AutomotiveEMS/electrical_load.mat', 'electrical_load', '-v7.3');
save('results/AutomotiveEMS/bus.mat', 'bus', '-v7.3');
save('results/AutomotiveEMS/uc_bank.mat', 'uc_bank', '-v7.3');
