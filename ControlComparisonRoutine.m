%% Inicializa modelo no Simulink

open_system('models/ControlComparison.slx', 'loadonly');

%% ParÃ¢metros temporais

T_s = 1e-6; % Passo de cÃ¡lculo utilizado pelo 'solver' local para sistemas fÃ­sicos [s]
T_k = 1e-4; % Passo de amostragem global de rotinas de controle [s]
t_f = 4e-0;	% Tempo total de simulaÃ§Ã£o [s]

%% Motor a combustÃ£o interna

% RazÃ£o entre a polia do motor a combustÃ£o e a polia do alternador
iceToAltRotRatio = 2.5;

% Pontos de interesse do perfil de velocidade
n_ice_i = 6e+3/iceToAltRotRatio;
n_ice_f = 2e+3/iceToAltRotRatio;
t_brake = 3e-0;                     % [s]
t_brake_i = (t_f - t_brake)/2;      % [s]
t_brake_f = t_brake_i + t_brake;    % [s]

%% Alternador

% Corrente de excitação máxima
i_f_max = 3.5;  % [A]

%% Retificador

% Filtro passivo
rectifier.filter.c = 10e-3;	% CapacitÃ¢ncia de filtro [F]

%% Carga elétrica

electrical_load.r = 1.0e+0;	% [Ohm]

%% Esquemas de controle

% Lista de títulos por esquema de controle
control_scheme_title = {'Esquema de controle convencional', ...
    'Esquema de controle baseado em excita{\c{c}}{\~{a}}o maximizada', ...
    'Esquema de controle baseado em rede neural'};

%
conventionalControlScheme = Simulink.Variant('control_scheme == 1');

%
fullExcitationCurrentControlScheme = Simulink.Variant('control_scheme == 2');

%
neuralControlScheme = Simulink.Variant('control_scheme == 3');
normalizationFactor = [4.0e+0 6.0e+3 2.0e+0];

%% ParÃ¢metros de simulaÃ§Ã£o

% ParÃ¢metros do 'solver' local para sistemas fÃ­sicos
set_param('ControlComparison/Solver Configuration', 'UseLocalSolver', 'on');
set_param('ControlComparison/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('ControlComparison/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('ControlComparison/Solver Configuration', 'DoFixedCost', 'on');
set_param('ControlComparison/Solver Configuration', 'MaxNonlinIter', '20');

% ParÃ¢metros do 'solver' global
set_param('ControlComparison', 'StopTime', num2str(t_f));

% Salva mundanÃ§as feitas no modelo
save_system('models/ControlComparison.slx');

%% ConfiguraÃ§Ã£o dos esquemas de controle como entrada do modelo no Simulink

% Inicialização da variável que define o esquema de controle selecionado
control_scheme = 1;

% 
for control_scheme_index = 1:length(control_scheme_title)
    simIn(control_scheme_index) = Simulink.SimulationInput('ControlComparison');
    simIn(control_scheme_index) = simIn(control_scheme_index).setVariable('control_scheme', control_scheme_index);
end

%% ExecuÃ§Ã£o da simulaÃ§Ã£o em ambiente Simulink

% ExecuÃ§Ã£o da simulaÃ§Ã£o paralelizada
simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
    'TransferBaseWorkspaceVariables', 'on');

%% Salva e finaliza modelo no Simulink

save_system('models/ControlComparison.slx');
close_system('models/ControlComparison.slx');

%% Registro de resultados obtidos na simulaÃ§Ã£o

for control_scheme_index = 1:length(control_scheme_title)
    % Motor a combustÃ£o interna
    ice.n{control_scheme_index} = simOut(control_scheme_index).n_ice;
    
    % Alternador
    alternator.rotor.n{control_scheme_index} = simOut(control_scheme_index).n_r;
    alternator.rotor.control.u{control_scheme_index} = simOut(control_scheme_index).u_i_f;
    alternator.rotor.l.i{control_scheme_index} = simOut(control_scheme_index).i_f;
    
    alternator.stator.input.e.value{control_scheme_index} = simOut(control_scheme_index).e_a_abc;
    alternator.stator.output.v{control_scheme_index} = simOut(control_scheme_index).v_a_abc;
    alternator.stator.output.i{control_scheme_index} = simOut(control_scheme_index).i_a_abc;
    
    % Retificador
    rectifier.control.u{control_scheme_index} = simOut(control_scheme_index).u_smr;
    
    % Carga
    electrical_load.v{control_scheme_index} = simOut(control_scheme_index).v_l;
    electrical_load.i{control_scheme_index} = simOut(control_scheme_index).i_l;
    electrical_load.z{control_scheme_index} = simOut(control_scheme_index).z_l;
    electrical_load.p{control_scheme_index} = simOut(control_scheme_index).p_l;
    
    % Bateria
    battery.v{control_scheme_index} = simOut(control_scheme_index).v_b;
end

%% Armazenamento dos resultados de simulaÃ§Ã£o

save('results/ControlComparison/ice.mat', 'ice', '-v7.3');
save('results/ControlComparison/alternator.mat', 'alternator', '-v7.3');
save('results/ControlComparison/rectifier.mat', 'rectifier', '-v7.3');
save('results/ControlComparison/electrical_load.mat', 'electrical_load', '-v7.3');
save('results/ControlComparison/battery.mat', 'battery', '-v7.3');

%% Montagem de séries temporais de valores esperados relativos aos MPPs

% Carregamento de mapa de dados de pontos de máxima potência
mpp_target = [];

try
    load('results/MPPTCurves/mpp_map.mat');
catch
    disp('MPP Map unavailable');
end

% Caso o arquivo tenha sido carregado com sucesso, executar a rotina
if (exist('mpp_u_3d', 'var'))
    % Inicialização dos pontos de interesse para comparação com as
    % estratégias de controle
    target_points = 20;
    t_target = round(linspace(t_brake_i, t_brake_f, target_points)', 5);
    i_f = zeros(target_points, 1);
    n_alt = zeros(target_points, 1);
    
    % Laço para determinar MPPs ideais para cada estratégia de controle
    for control_scheme_index = 1:length(control_scheme_title)
        % Armazenamento das variáveis de interesse
        t = round(alternator.rotor.l.i{control_scheme_index}.time, 10);
        i_f = alternator.rotor.l.i{control_scheme_index}.data;
        n_alt = alternator.rotor.n{control_scheme_index}.data;
        
        % Determinação dos índices dos pontos de interesse
        target_elements = ismember(t, t_target);
        
        t = t(target_elements);
        i_f = i_f(target_elements);
        n_alt = n_alt(target_elements);
        
        % Remoção de pontos idênticos
        [~, target_elements, ~] = unique(t, 'first');
        
        t = t(target_elements);
        i_f = i_f(target_elements);
        n_alt = n_alt(target_elements);
        r_l = electrical_load.r.*ones(target_points, 1);
        
        % Determinação dos MPPs ideais via interpolação
        mpp_target_u_tmp = interpn(i_f_list, n_alt_list, r_l_list, mpp_u_3d, i_f, n_alt, r_l, 'spline', -1);
        mpp_target_p_tmp = interpn(i_f_list, n_alt_list, r_l_list, mpp_p_3d, i_f, n_alt, r_l, 'spline', -1);
        
        % Formação das séries temporais dos MPPs ideais
        mpp_target_u_tmp = timeseries(mpp_target_u_tmp, t);
        mpp_target_p_tmp = timeseries(mpp_target_p_tmp, t);
        
        % Armazenamento das variáveis para traço de figuras
        mpp_target_u{control_scheme_index} = mpp_target_u_tmp;
        mpp_target_p{control_scheme_index} = mpp_target_p_tmp;
    end
end

%% Análise de potência e ciclo de trabalho para cada esquema de controle

% Índice de figuras
figure_index = 0;

% Caso o arquivo tenha sido carregado com sucesso, executar a rotina
if (exist('mpp_u_3d', 'var'))
    
    % Laço de traço de figuras
    for control_scheme_index = 1:length(control_scheme_title)
        figure_index = figure_index + 1;
        figure(figure_index)
        
        subplot(4, 1, 1)
        plot(electrical_load.p{control_scheme_index}, 'b-', 'DisplayName', '$P_{l}$');
        hold on;
        plot(mpp_target_p{control_scheme_index}, 'ro', 'DisplayName', '$P_{MPP}^{*}$');
        hold off;
        title('Pot{\^{e}}ncia el{\''{e}}trica consumida pela carga');
        xlabel('$t$ [s]');
        ylabel('$P_{l}$ [W]');
        legend('show');
        grid on;
        
        subplot(4, 1, 2)
        plot(rectifier.control.u{control_scheme_index}, 'b-', 'DisplayName', '$d_{SMR}$');
        hold on;
        plot(mpp_target_u{control_scheme_index}, 'ro', 'DisplayName', '$d_{MPP}^{*}$');
        hold off;
        title('Ciclo de trabalho das chaves do retificador semi-controlado');
        xlabel('$t$ [s]');
        ylabel('$d_{SMR}$');
        legend('show');
        grid on;
        
        subplot(4, 1, 3)
        plot(alternator.rotor.l.i{control_scheme_index});
        title('Corrente de excita{\c{c}}{\~{a}}o do alternador');
        xlabel('$t$ [s]');
        ylabel('$i_{f} [A]$');
        grid on;
        
        subplot(4, 1, 4)
        plot(alternator.rotor.n{control_scheme_index});
        title('Velocidade do alternador');
        xlabel('$t$ [s]');
        ylabel('$n_{alt} [rpm]$');
        grid on;
        
        suptitle(control_scheme_title{control_scheme_index});
    end
end

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = sprintf('results/ControlComparison/Figura_%d', i);
    saveFigure(figure(i), fileName, 'fig');
end
