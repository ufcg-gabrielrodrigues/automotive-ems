%% Inicializa modelo no Simulink

open_system('models/ThermalAnalysis.slx', 'loadonly');

%% Parâmetros temporais

T_s = 1e-6;     % Passo de cálculo utilizado pelo 'solver' [s]
t_f = 1.0e-2;   % Tempo total de simulação [s]

%% Varredura de parâmetros

% Lista de parâmetros a serem varridos individualmente
i_f_list = [1.0 3.0 5.0]';                  % Corrente de excitação máxima [A]
n_r_list = [2000 3500 5000 7500]';          % Velocidade do alternador [rpm]
T_list = [20.0 50.0 100.0 150.0 200.0]';    % Temperaturas da resistência de enrolamento de estator [oC]
v_o_list = (0.0:1.0:80.0)';                 % Tensão de saída [V]

%% Parâmetros auxiliares para figuras

% Índice de figuras
figure_index = 0;

% Cores
colors = distinguishable_colors(length(T_list));

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
set_param('ThermalAnalysis/Solver Configuration', 'UseLocalSolver', 'on');
set_param('ThermalAnalysis/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('ThermalAnalysis/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('ThermalAnalysis/Solver Configuration', 'DoFixedCost', 'on');
set_param('ThermalAnalysis/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
set_param('ThermalAnalysis', 'StopTime', num2str(t_f));

% Salva mundanças feitas no modelo
save_system('models/ThermalAnalysis.slx');

%% Configuração dos casos de teste para carga de tensão constante

% Configuração dos casos de teste como grid 4D
[i_f_grid, n_r_grid, T_grid, v_o_grid] = ndgrid(i_f_list, n_r_list, T_list, v_o_list);

% Configuração dos casos de teste como entrada do modelo no Simulink

for i_f_index = 1:length(i_f_list)
    i_f = i_f_grid(i_f_index, 1, 1, 1);
    
    for n_r_index = 1:length(n_r_list)
        n_r = n_r_grid(1, n_r_index, 1, 1);
        
        for T_index = 1:length(T_list)
            T = T_grid(1, 1, T_index, 1);
            
            for v_o_index = 1:length(v_o_list)
                v_o = v_o_grid(1, 1, 1, v_o_index);
                
                simIn(i_f_index, n_r_index, T_index, v_o_index) = Simulink.SimulationInput('ThermalAnalysis');
                
                if (alternator.stator.connection == y)
                    simIn(i_f_index, n_r_index, T_index, v_o_index) = simIn(i_f_index, n_r_index, T_index, v_o_index).setBlockParameter('ThermalAnalysis/Alternator/Armature circuit/Armature circuit (Y)/r_a', 'T', num2str(T));
                    simIn(i_f_index, n_r_index, T_index, v_o_index) = simIn(i_f_index, n_r_index, T_index, v_o_index).setBlockParameter('ThermalAnalysis/Alternator/Armature circuit/Armature circuit (Y)/r_b', 'T', num2str(T));
                    simIn(i_f_index, n_r_index, T_index, v_o_index) = simIn(i_f_index, n_r_index, T_index, v_o_index).setBlockParameter('ThermalAnalysis/Alternator/Armature circuit/Armature circuit (Y)/r_c', 'T', num2str(T));
                elseif (alternator.stator.connection == delta)
                    simIn(i_f_index, n_r_index, T_index, v_o_index) = simIn(i_f_index, n_r_index, T_index, v_o_index).setBlockParameter('ThermalAnalysis/Alternator/Armature circuit/Armature circuit (Delta)/r_a', 'T', num2str(T));
                    simIn(i_f_index, n_r_index, T_index, v_o_index) = simIn(i_f_index, n_r_index, T_index, v_o_index).setBlockParameter('ThermalAnalysis/Alternator/Armature circuit/Armature circuit (Delta)/r_b', 'T', num2str(T));
                    simIn(i_f_index, n_r_index, T_index, v_o_index) = simIn(i_f_index, n_r_index, T_index, v_o_index).setBlockParameter('ThermalAnalysis/Alternator/Armature circuit/Armature circuit (Delta)/r_c', 'T', num2str(T));
                end
                
                simIn(i_f_index, n_r_index, T_index, v_o_index) = simIn(i_f_index, n_r_index, T_index, v_o_index).setBlockParameter('ThermalAnalysis/i_f', 'Value', num2str(i_f));
                simIn(i_f_index, n_r_index, T_index, v_o_index) = simIn(i_f_index, n_r_index, T_index, v_o_index).setBlockParameter('ThermalAnalysis/n_r', 'Value', num2str(n_r));
                simIn(i_f_index, n_r_index, T_index, v_o_index) = simIn(i_f_index, n_r_index, T_index, v_o_index).setBlockParameter('ThermalAnalysis/v_o', 'DC', num2str(v_o));
            end
        end
    end
end

% Transformação de matriz de entradas em vetor
simIn = reshape(simIn, [length(i_f_list)*length(n_r_list)*length(T_list)*length(v_o_list) 1]);

%% Análise do efeito da variação da tensão na carga

% Execução da simulação paralelizada
simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
    'TransferBaseWorkspaceVariables', 'on');

simIn = reshape(simIn, [length(i_f_list) length(n_r_list) length(T_list) length(v_o_list)]);
simOut = reshape(simOut, [length(i_f_list) length(n_r_list) length(T_list) length(v_o_list)]);

for i_f_index = 1:length(i_f_list)
    for n_r_index = 1:length(n_r_list)
        for T_index = 1:length(T_list)
            for v_o_index = 1:length(v_o_list)
                if (isempty(simOut(i_f_index, n_r_index, T_index, v_o_index).ErrorMessage))
                    p_o(i_f_index, n_r_index, T_index, v_o_index) = mean(simOut(i_f_index, n_r_index, T_index, v_o_index).p_l.data(round(end/2):end));
                else
                    p_o(i_f_index, n_r_index, T_index, v_o_index) = nan;
                end
            end
        end
    end
end

%% Identificação dos pontos de máxima potência indexados pela tensão da carga

% 
[p_o_mpp, v_o_mpp_index] = max(p_o, [], 4);
v_o_mpp = v_o_list(v_o_mpp_index);

%% Traço dos resultados relativos à variação da tensão da carga

% 
tit = cell(length(i_f_list) * length(n_r_list));

% 
for i_f_index = 1:length(i_f_list)
    
    for n_r_index = 1:length(n_r_list)
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for T_index = 1:length(T_list)
            
            plot(v_o_list, squeeze(p_o(i_f_index, n_r_index, T_index, :)), '-', 'Color', colors(T_index, :), ...
                'DisplayName', ['$T = ' num2str(T_list(T_index)) '^{o}\textrm{C}$']);
            hold on;
            plot(v_o_mpp(i_f_index, n_r_index, T_index, :), p_o_mpp(i_f_index, n_r_index, T_index, :), 'o', ...
                'Color', colors(T_index, :), 'HandleVisibility', 'off');
            
            legend('off');
            legend('show');
            
        end
        
        tit{figure_index} = ['po-vo-if-' num2str(i_f_list(i_f_index)) '-nr-' num2str(n_r_list(n_r_index))];
        xlabel('$v_{o}\,[\textrm{V}]$');
        ylabel('$p_{o}\,[\textrm{W}]$');
        leg = legend;
        leg.Location = 'NorthWest';
        grid on;
    end
end

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = ['results/ThermalAnalysis/' tit{i}];
    saveFigure(figure(i), fileName, 'fig');
    saveFigure(figure(i), fileName, 'eps');
end

%% Armazenamento dos resultados de simulação

save('results/ThermalAnalysis/p_o.mat', 'simIn', 'simOut', 'p_o', 'p_o_mpp', 'v_o_mpp', '-v7.3');
save('results/ThermalAnalysis/simEnv.mat', 'alternator', 'rectifier', ...
    'i_f_list', 'n_r_list', 'T_list', 'v_o_list', '-v7.3');
