%% Parâmetros temporais

T_s = 1e-6; % Passo de cálculo utilizado pelo 'solver' local para sistemas físicos [s]
T_k = 1e-4; % Passo de amostragem global de rotinas de controle [s]
t_f = 5.0;  % Tempo total de simulação [s]

%% Motor a combustão interna

% Razão entre a polia do motor a combustão e a polia do alternador
iceToAltRotRatio = 2.5;

% Pontos de interesse do perfil de velocidade
n_ice_i = 6e+3/iceToAltRotRatio;
n_ice_f = 2e+3/iceToAltRotRatio;
t_brake_i = t_f*0.20;
t_brake_f = t_f*0.80;

%% Alternador

% Corrente de excita��o m�xima
i_f_max = 3.0;  % [A]

%% Retificador

% Filtro passivo
rectifier.filter.c = 10e-3;	% Capacitância de filtro [F]

%% Carga el�trica

electrical_load.r = 0.5;

%% Esquemas de controle

% Lista de t�tulos por esquema de controle
control_scheme_title = {'Esquema de controle convencional', ...
    'Esquema de controle baseado em excita{\c{c}}{\~{a}}o maximizada', ...
    'Esquema de controle baseado em rede neural', ...
    'Esquema de controle baseado em superf{\''{i}}cie ajustada'};

%
conventionalControlScheme = Simulink.Variant('control_scheme == 1');

%
fullExcitationCurrentControlScheme = Simulink.Variant('control_scheme == 2');

%
neuralControlScheme = Simulink.Variant('control_scheme == 3');

%
fittedControlScheme = Simulink.Variant('control_scheme == 4');
load('src/control_scheme/fitted_mpp_surface/fittedMPPSurface.mat');

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

for control_scheme_index = 1:length(control_scheme_title)
    fprintf('Running control scheme #%d...\n', control_scheme_index);
    control_scheme = control_scheme_index;
    simout{control_scheme_index} = sim('ControlComparison', simulationParameters);
end

%% Salva e finaliza modelo no Simulink

save_system('models/ControlComparison.slx');
close_system('models/ControlComparison.slx');

%% Registro de resultados obtidos na simulação

for control_scheme_index = 1:length(control_scheme_title)
    % Motor a combustão interna
    ice.n{control_scheme_index} = simout{control_scheme_index}.n_ice;
    
    % Alternador
    alternator.rotor.n{control_scheme_index} = simout{control_scheme_index}.n_r;
    alternator.rotor.control.u{control_scheme_index} = simout{control_scheme_index}.u_i_f;
    alternator.rotor.l.i{control_scheme_index} = simout{control_scheme_index}.i_f;
    
    alternator.stator.input.e.value{control_scheme_index} = simout{control_scheme_index}.e_a_abc;
    alternator.stator.output.v{control_scheme_index} = simout{control_scheme_index}.v_a_abc;
    alternator.stator.output.i{control_scheme_index} = simout{control_scheme_index}.i_a_abc;
    
    % Retificador
    rectifier.control.u{control_scheme_index} = simout{control_scheme_index}.u_smr;
    
    % Carga
    electrical_load.v{control_scheme_index} = simout{control_scheme_index}.v_l;
    electrical_load.i{control_scheme_index} = simout{control_scheme_index}.i_l;
    electrical_load.z{control_scheme_index} = simout{control_scheme_index}.z_l;
    electrical_load.p{control_scheme_index} = simout{control_scheme_index}.p_l;
    
    % Bateria
    battery.v{control_scheme_index} = simout{control_scheme_index}.v_b;
end

%% Armazenamento dos resultados de simulação

save('results/ControlComparison/ice.mat', 'ice', '-v7.3');
save('results/ControlComparison/alternator.mat', 'alternator', '-v7.3');
save('results/ControlComparison/rectifier.mat', 'rectifier', '-v7.3');
save('results/ControlComparison/electrical_load.mat', 'electrical_load', '-v7.3');
save('results/ControlComparison/battery.mat', 'battery', '-v7.3');

%% Montagem de s�ries temporais de valores esperados relativos aos MPPs

% Carregamento de matrix de dados de pontos de m�xima pot�ncia
mpp_target = [];

try
    load('results/MPPTCurves/mpp_matrix.mat');
    mpp_target = mpp_matrix(mpp_matrix(:, 3) == electrical_load.r, :);
catch
    disp('MPP Matrix unavailable');
end

% Caso a matriz tenha sido carregada com sucesso, executar a rotina
if (~isempty(mpp_target))
    % Remo��o de valores fora da regi�o de interesse de velocidade
    n_r_min = min(iceToAltRotRatio*[n_ice_i n_ice_f]);
    n_r_max = max(iceToAltRotRatio*[n_ice_i n_ice_f]);
    
    mpp_target((mpp_target(:, 2) < n_r_min) | (mpp_target(:, 2) > n_r_max), :) = [];
    
    % Determina��o da equa��o da reta que descreve a varia��o de velocidade
    a = iceToAltRotRatio*(n_ice_f - n_ice_i)/(t_brake_f - t_brake_i);
    b = iceToAltRotRatio*n_ice_i - a*t_brake_i;
    
    % Determina��o dos instantes correspondente aos valores de interesse da
    % velocidade
    time = (mpp_target(:, 2) - b)/a;
    
    % Ordena��o da matriz de casos esperados com rela��o ao tempo
    mpp_target = [time mpp_target];
    [~, sorted_indexes] = sort(mpp_target(:, 1), 'ascend');
    mpp_target = mpp_target(sorted_indexes, :);
    
    % Montagem das s�ries temporais relativas ao ciclo de trabalho e �
    % pot�ncia
    mpp_target_u = timeseries(mpp_target(:, 5), mpp_target(:, 1));
    mpp_target_p = timeseries(mpp_target(:, 6), mpp_target(:, 1));
end

%% An�lise de pot�ncia e ciclo de trabalho para cada esquema de controle

% �ndice de figuras
figure_index = 0;
    
% La�o de tra�o de figuras
for control_scheme_index = 1:length(control_scheme_title)
    figure_index = figure_index + 1;
    figure(figure_index)
    
    subplot(2, 1, 1)
    plot(electrical_load.p{control_scheme_index}, 'b-', 'DisplayName', '$P_{l}$');
    hold on;
    plot(mpp_target_p, 'ro', 'DisplayName', '$P_{MPP}^{*}$');
    hold off;
    title('Pot{\^{e}}ncia el{\''{e}}trica consumida pela carga');
    xlabel('$t$ [s]');
    ylabel('$P_{l}$ [W]');
    legend('show');
    grid on;
    
    subplot(2, 1, 2)
    plot(rectifier.control.u{control_scheme_index}, 'b-', 'DisplayName', '$d_{SMR}$');
    hold on;
    plot(mpp_target_u, 'ro', 'DisplayName', '$d_{MPP}^{*}$');
    hold off;
    title('Ciclo de trabalho das chaves do retificador semi-controlado');
    xlabel('$t$ [s]');
    ylabel('$d_{SMR}$');
    legend('show');
    grid on;
    
    suptitle(control_scheme_title{control_scheme_index});
end

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = sprintf('results/ControlComparison/Figura_%d', i);
    saveFigure(figure(i), fileName, 'fig');
end
