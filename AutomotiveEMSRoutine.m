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

% Parâmetros de esquema de controle híbrido
hybrid_control_cases = [1.0 0.0; 0.0 10.0; 1.0 0.1];

%% Parâmetros auxiliares para figuras

[num_cases, ~] = size(hybrid_control_cases);

% �?ndice de figuras
figure_index = 0;

% Cores
colors = distinguishable_colors(num_cases + 1);

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

%% Conversor Buck

% Eficiência
buck.efficiency = 0.85;
buck.v_o_max = 13.5;

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

%% Configuração dos casos de teste como entrada do modelo no Simulink

for test_case_index = 1:num_cases
    test_case = hybrid_control_cases(test_case_index, :);
    
    hybrid_control.k_p(test_case_index) = test_case(1);
    hybrid_control.k_i(test_case_index) = test_case(2);
    
    simIn(test_case_index) = Simulink.SimulationInput('AutomotiveEMS');
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('AutomotiveEMS/Hybrid Load Matching Controller Scheme/Hybrid Load Matching Controller/k_p', 'Gain', num2str(hybrid_control.k_p(test_case_index)));
    simIn(test_case_index) = simIn(test_case_index).setBlockParameter('AutomotiveEMS/Hybrid Load Matching Controller Scheme/Hybrid Load Matching Controller/k_i', 'Gain', num2str(hybrid_control.k_i(test_case_index)));
end

%% Execução da simulação em ambiente Simulink

% Execução da simulação paralelizada
simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
    'TransferBaseWorkspaceVariables', 'on');

%% Redefinição de parâmetros de alternador

% Atualização de parâmetro para valor padrão: fator de acoplamento
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', 'AutomotiveEMS/Hybrid Load Matching Controller Scheme/Hybrid Load Matching Controller/Load Matching Switched-Mode Rectifier Controller/MATLAB Function');
blockHandle.Script = strrep(blockHandle.Script, k_e_local, k_e_default_local);

%% Salva e finaliza modelo no Simulink

save_system('models/AutomotiveEMS.slx');
close_system('models/AutomotiveEMS.slx');

%% Registro de resultados obtidos no caso de teste

% Laço de iterações por casos de teste
for test_case_index = 1:num_cases
    % Motor a combustão interna
    ice.n(test_case_index) = simOut(test_case_index).n_ice;
    
    % Alternador
    alternator.rotor.control.u(test_case_index) = simOut(test_case_index).u_i_f;
    alternator.rotor.l.i(test_case_index) = simOut(test_case_index).i_f;
    
    alternator.stator.input.e(test_case_index) = simOut(test_case_index).e_a_abc;
    alternator.stator.output.v(test_case_index) = simOut(test_case_index).v_a_abc;
    alternator.stator.output.i(test_case_index) = simOut(test_case_index).i_a_abc;
    
    % Retificador
    rectifier.control.u_lm(test_case_index) = simOut(test_case_index).u_lm;
    rectifier.control.u_esc(test_case_index) = simOut(test_case_index).u_esc;
    rectifier.control.u_smr(test_case_index) = simOut(test_case_index).u_smr;
    
    % Carga elétrica
    electrical_load.v(test_case_index) = simOut(test_case_index).v_l;
    electrical_load.i(test_case_index) = simOut(test_case_index).i_l;
    electrical_load.p(test_case_index) = simOut(test_case_index).p_l;
    
    % Barramento de circuito secundário
    bus.v(test_case_index) = simOut(test_case_index).v_bus;
    bus.i(test_case_index) = simOut(test_case_index).i_bus;
    bus.p(test_case_index) = simOut(test_case_index).p_bus;
    
    % Banco de supercapacitores
    uc_bank.v(test_case_index) = simOut(test_case_index).v_uc;
    uc_bank.i(test_case_index) = simOut(test_case_index).i_uc;
    uc_bank.p(test_case_index) = simOut(test_case_index).p_uc;
end

%% 

% Carregamento de superfície ajustada de máxima potência
try
    load('results/PowerAnalysis/P_v_o.mat', 'P_v_o_mpp_sim_fit');
catch
    disp('Maximum power fit unavailable');
end

%% Análise de potência e ciclo de trabalho para cada esquema de controle

% 
figure_index = figure_index + 1;
figure(figure_index)

for test_case_index = 1:num_cases
    plot(uc_bank.p(test_case_index).Time, uc_bank.p(test_case_index).Data, 'Color', colors(test_case_index, :), ...
        'DisplayName', ['$k_{p} = ' num2str(hybrid_control.k_p(test_case_index)) '$; $k_{i} = ' num2str(hybrid_control.k_i(test_case_index)) '$']);
    hold on;
end

if (exist('P_v_o_mpp_sim_fit', 'var'))
    P_v_o_mpp = buck.efficiency*P_v_o_mpp_sim_fit(iceToAltRotRatio*ice.n(test_case_index).Data, alternator.rotor.l.i(test_case_index).Data);
    P_v_o_mpp = timeseries(P_v_o_mpp, ice.n(test_case_index).Time);
    
    plot(P_v_o_mpp.Time, P_v_o_mpp.Data, 'Color', colors(num_cases + 1, :), ...
        'DisplayName', 'Pot{\^{e}}ncia m{\''{a}}xima');
end

hold off;
legend('off');
legend('show');

title('Pot{\^{e}}ncia transferida ao banco de supercapacitores');
xlabel('$t$ $[s]$');
ylabel('$P_{uc}$ $[W]$');
leg = legend;
leg.Location = 'NorthEast';
grid on;

% 
figure_index = figure_index + 1;
figure(figure_index)

for test_case_index = 1:num_cases
    plot(ice.n(test_case_index).Time, ice.n(test_case_index).Data*iceToAltRotRatio, 'Color', colors(test_case_index, :), ...
        'DisplayName', ['$k_{p} = ' num2str(hybrid_control.k_p(test_case_index)) '$; $k_{i} = ' num2str(hybrid_control.k_i(test_case_index)) '$']);
    hold on;
end

hold off;
legend('off');
legend('show');

title('Velocidade do alternador');
xlabel('$t$ $[s]$');
ylabel('$n_{r}$ $[rpm]$');
leg = legend;
leg.Location = 'NorthEast';
grid on;

% 
figure_index = figure_index + 1;
figure(figure_index)

for test_case_index = 1:num_cases
    plot(alternator.rotor.l.i(test_case_index).Time, alternator.rotor.l.i(test_case_index).Data, 'Color', colors(test_case_index, :), ...
        'DisplayName', ['$k_{p} = ' num2str(hybrid_control.k_p(test_case_index)) '$; $k_{i} = ' num2str(hybrid_control.k_i(test_case_index)) '$']);
    hold on;
end

hold off;
legend('off');
legend('show');

title('Corrente de excita{\c{c}}{\~{a}}o do alternador');
xlabel('$t$ $[s]$');
ylabel('$i_{f}$ $[A]$');
leg = legend;
leg.Location = 'NorthEast';
grid on;

% 
figure_index = figure_index + 1;
figure(figure_index)

for test_case_index = 1:num_cases
    plot(bus.v(test_case_index).Time, bus.v(test_case_index).Data, 'Color', colors(test_case_index, :), ...
        'DisplayName', ['$k_{p} = ' num2str(hybrid_control.k_p(test_case_index)) '$; $k_{i} = ' num2str(hybrid_control.k_i(test_case_index)) '$']);
    hold on;
end

hold off;
legend('off');
legend('show');

title('Tens{\~{a}}o de barramento do circuito secund{\''{a}}rio');
xlabel('$t$ $[s]$');
ylabel('$v_{dc}$ $[V]$');
leg = legend;
leg.Location = 'NorthEast';
grid on;

% 
figure_index = figure_index + 1;
figure(figure_index)

subplot(3, 1, 1)
for test_case_index = 1:num_cases
    plot(rectifier.control.u_lm(test_case_index).Time, rectifier.control.u_lm(test_case_index).Data, 'Color', colors(test_case_index, :), ...
        'DisplayName', ['$k_{p} = ' num2str(hybrid_control.k_p(test_case_index)) '$; $k_{i} = ' num2str(hybrid_control.k_i(test_case_index)) '$']);
    hold on;
end

hold off;
legend('off');
legend('show');

title('Sinal de controle baseado em modelo');
xlabel('$t$ $[s]$');
ylabel('$u$');
leg = legend;
leg.Location = 'NorthEast';
grid on;

subplot(3, 1, 2)
for test_case_index = 1:num_cases
    plot(rectifier.control.u_esc(test_case_index).Time, rectifier.control.u_esc(test_case_index).Data, 'Color', colors(test_case_index, :), ...
        'DisplayName', ['$k_{p} = ' num2str(hybrid_control.k_p(test_case_index)) '$; $k_{i} = ' num2str(hybrid_control.k_i(test_case_index)) '$']);
    hold on;
end

hold off;
legend('off');
legend('show');

title('Sinal de controle baseado em otimiza{\c{c}}{\~{a}}o');
xlabel('$t$ $[s]$');
ylabel('$u$');
leg = legend;
leg.Location = 'NorthEast';
grid on;

subplot(3, 1, 3)
for test_case_index = 1:num_cases
    plot(rectifier.control.u_smr(test_case_index).Time, rectifier.control.u_smr(test_case_index).Data, 'Color', colors(test_case_index, :), ...
        'DisplayName', ['$k_{p} = ' num2str(hybrid_control.k_p(test_case_index)) '$; $k_{i} = ' num2str(hybrid_control.k_i(test_case_index)) '$']);
    hold on;
end

hold off;
legend('off');
legend('show');

title('Sinal de controle resultante');
xlabel('$t$ $[s]$');
ylabel('$u$');
leg = legend;
leg.Location = 'NorthEast';
grid on;

suptitle('Sinais de controle de retificador');

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = sprintf('results/AutomotiveEMS/Figura_%d', i);
    saveFigure(figure(i), fileName, 'fig');
end

%% Armazenamento dos resultados de simulação

save('results/AutomotiveEMS/ice.mat', 'ice', '-v7.3');
save('results/AutomotiveEMS/alternator.mat', 'alternator', '-v7.3');
save('results/AutomotiveEMS/rectifier.mat', 'rectifier', '-v7.3');
save('results/AutomotiveEMS/battery.mat', 'battery', '-v7.3');
save('results/AutomotiveEMS/electrical_load.mat', 'electrical_load', '-v7.3');
save('results/AutomotiveEMS/bus.mat', 'bus', '-v7.3');
save('results/AutomotiveEMS/uc_bank.mat', 'uc_bank', '-v7.3');
save('results/AutomotiveEMS/hybrid_control.mat', 'hybrid_control', '-v7.3');
save('results/AutomotiveEMS/buck.mat', 'buck', '-v7.3');
