%% Inicializa modelo no Simulink

open_system('models/PowerAnalysis.slx', 'loadonly');

%% Parâmetros temporais

T_s = 1e-6;     % Passo de cálculo utilizado pelo 'solver' [s]
t_f = 1.0e-2;   % Tempo total de simulação [s]

%% Modelos de carga

% Lista de títulos por modelo de carga
load_model_title = {'Carga de tens{\~{a}}o constante', 'Carga de imped{\^{a}}ncia constante'};

% Carga de tensão constante
constantVoltageLoad = Simulink.Variant('load_model == 1');

% Carga de impedância constante
constantImpedanceLoad = Simulink.Variant('load_model == 2');

% Inicialização da variável de escolha do modelo
load_model = 1;

%% Alternador

% Efeito térmico na resistência do circuito de estator
T = 150;	% [oC]
alternator.stator.r.value = alternator.stator.r.function(T);

% Fator de acoplamento
if (isfield(alternator.k_e, 'function'))
    k_e_fun = @(i_f) alternator.k_e.function(i_f);
else
    k_e_fun = @(i_f) alternator.k_e.value;
end

if (alternator.stator.connection == delta)
    k_e_fun = @(i_f) k_e_fun(i_f)./sqrt(3);
end

% Indutância de estator
if (isfield(alternator.stator.l, 'function'))
    l_s_fun = @(i_f) alternator.stator.l.function(i_f);
else
    l_s_fun = @(i_f) alternator.stator.l.value;
end

if (alternator.stator.connection == delta)
    l_s_fun = @(i_f) l_s_fun(i_f)./3;
end

% Função de cálculo da frequência elétrica
omega_e = @(n_r) n_r.*(2.*pi./60).*alternator.p;

% Função de cálculo da tensão induzida no estator
v_s = @(n_r, i_f) k_e_fun(i_f).*omega_e(n_r).*i_f;

%% Modelos analíticos para cálculo de potência

% Carga de tensão constante
P_v_o = @(n_r, i_f, v_o) (3.*v_o./pi).*(sqrt(v_s(n_r, i_f).^2 - (2.*v_o./pi).^2))./(omega_e(n_r).*l_s_fun(i_f));

% Carga de impedância constante
P_z_o = @(n_r, i_f, z_o) ((3.*pi.*v_s(n_r, i_f)).^2.*z_o)./((pi.^2.*omega_e(n_r).*l_s_fun(i_f)).^2 + (6.*z_o).^2);

%% Varredura de parâmetros

% Lista de parâmetros a serem varridos individualmente
i_f_list = [0.01 0.5:0.5:5.0]';      % Corrente de excitação máxima [A]
n_r_list = (2000:500:7500)';       	% Velocidade do alternador [rpm]
v_o_list = (0.0:1.0:80.0)';       	% Tensão de saída [V]
z_o_list = [0.01 0.05:0.05:2.0]';   % Impedância de saída [Ohm]

%% Parâmetros auxiliares para figuras

% Índice de figuras
figure_index = 0;

% Cores
colors_n_r = distinguishable_colors(length(n_r_list));
colors_i_f = distinguishable_colors(length(i_f_list));

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
set_param('PowerAnalysis/Solver Configuration', 'UseLocalSolver', 'on');
set_param('PowerAnalysis/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('PowerAnalysis/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('PowerAnalysis/Solver Configuration', 'DoFixedCost', 'on');
set_param('PowerAnalysis/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
set_param('PowerAnalysis', 'StopTime', num2str(t_f));

% Salva mundanças feitas no modelo
save_system('models/PowerAnalysis.slx');

%% Configuração dos casos de teste para carga de tensão constante

% Configuração dos casos de teste como entrada do modelo analítico
[i_f_grid, n_r_grid, v_o_grid] = meshgrid(i_f_list, n_r_list, v_o_list);

% Configuração dos casos de teste como entrada do modelo no Simulink
load_model = 1;
clear simIn;

for i_f_index = 1:length(i_f_list)
    i_f = i_f_grid(1, i_f_index, 1);
    
    for n_r_index = 1:length(n_r_list)
        n_r = n_r_grid(n_r_index, 1, 1);
        
        for v_o_index = 1:length(v_o_list)
            v_o = v_o_grid(1, 1, v_o_index);
            
            simIn(n_r_index, i_f_index, v_o_index) = Simulink.SimulationInput('PowerAnalysis');
            simIn(n_r_index, i_f_index, v_o_index) = simIn(n_r_index, i_f_index, v_o_index).setVariable('load_model', load_model);
            simIn(n_r_index, i_f_index, v_o_index) = simIn(n_r_index, i_f_index, v_o_index).setBlockParameter('PowerAnalysis/i_f', 'Value', num2str(i_f));
            simIn(n_r_index, i_f_index, v_o_index) = simIn(n_r_index, i_f_index, v_o_index).setBlockParameter('PowerAnalysis/n_r', 'Value', num2str(n_r));
            simIn(n_r_index, i_f_index, v_o_index) = simIn(n_r_index, i_f_index, v_o_index).setBlockParameter('PowerAnalysis/Load/Voltage/v_o', 'DC', num2str(v_o));
        end
    end
end

% Transformação de matriz de entradas em vetor
simIn = reshape(simIn, [length(n_r_list)*length(i_f_list)*length(v_o_list) 1]);

%% Análise do efeito da variação da tensão na carga

% Modelo analítico
P_v_o_ana = P_v_o(n_r_grid, i_f_grid, v_o_grid);
P_v_o_ana(imag(P_v_o_ana) ~= 0) = 0;

% Execução da simulação paralelizada
simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
    'TransferBaseWorkspaceVariables', 'on');

simIn = reshape(simIn, [length(n_r_list) length(i_f_list) length(v_o_list)]);
simOut = reshape(simOut, [length(n_r_list) length(i_f_list) length(v_o_list)]);

for i_f_index = 1:length(i_f_list)
    for n_r_index = 1:length(n_r_list)
        for v_o_index = 1:length(v_o_list)
            P_v_o_sim(n_r_index, i_f_index, v_o_index) = mean(simOut(n_r_index, i_f_index, v_o_index).p_l.data(round(end/2):end));
        end
    end
end

%% Identificação dos pontos de máxima potência indexados pela tensão da carga

% 
[P_v_o_mpp_ana, v_o_mpp_index_ana] = max(P_v_o_ana, [], 3);
v_o_mpp_ana = v_o_list(v_o_mpp_index_ana);

% 
[P_v_o_mpp_sim, v_o_mpp_index_sim] = max(P_v_o_sim, [], 3);
v_o_mpp_sim = v_o_list(v_o_mpp_index_sim);

%% Traço dos resultados relativos à variação da tensão da carga

% 
for i_f_index = 1:length(i_f_list)
    
    figure_index = figure_index + 1;
    figure(figure_index)
    
    for n_r_index = 1:length(n_r_list)
        
        plot(v_o_list, squeeze(P_v_o_ana(n_r_index, i_f_index, :)), '-', 'Color', colors_n_r(n_r_index, :), ...
            'DisplayName', ['$n_{r} = ' num2str(n_r_list(n_r_index)) ' rpm$']);
        hold on;
        plot(v_o_mpp_ana(n_r_index, i_f_index, :), P_v_o_mpp_ana(n_r_index, i_f_index, :), 'o', ...
            'Color', colors_n_r(n_r_index, :), 'HandleVisibility', 'off');
        hold on;
        plot(v_o_list, squeeze(P_v_o_sim(n_r_index, i_f_index, :)), '--', 'Color', colors_n_r(n_r_index, :), ...
            'HandleVisibility', 'off');
        hold on;
        plot(v_o_mpp_sim(n_r_index, i_f_index, :), P_v_o_mpp_sim(n_r_index, i_f_index, :), 'o', ...
            'Color', colors_n_r(n_r_index, :), 'HandleVisibility', 'off');
        
        legend('off');
        legend('show');
    end
    
    title(['Curvas de pot{\^{e}}ncia indexadas pela tens{\~{a}}o da carga ($i_{f}$ $=$ $' num2str(i_f_list(i_f_index)) '$ $[A]$)']);
    xlabel('$v_o$ $[V]$');
    ylabel('$P_o$ $[W]$');
    leg = legend;
    leg.Location = 'NorthWest';
    grid on;
end

% 
for n_r_index = 1:length(n_r_list)
    
    figure_index = figure_index + 1;
    figure(figure_index)
    
    for i_f_index = 1:length(i_f_list)
        
        plot(v_o_list, squeeze(P_v_o_ana(n_r_index, i_f_index, :)), '-', 'Color', colors_i_f(i_f_index, :), ...
            'DisplayName', ['$i_{f} = ' num2str(i_f_list(i_f_index)) ' A$']);
        hold on;
        plot(v_o_mpp_ana(n_r_index, i_f_index, :), P_v_o_mpp_ana(n_r_index, i_f_index, :), 'o', ...
            'Color', colors_i_f(i_f_index, :), 'HandleVisibility', 'off');
        hold on;
        plot(v_o_list, squeeze(P_v_o_sim(n_r_index, i_f_index, :)), '--', 'Color', colors_i_f(i_f_index, :), ...
            'HandleVisibility', 'off');
        hold on;
        plot(v_o_mpp_sim(n_r_index, i_f_index, :), P_v_o_mpp_sim(n_r_index, i_f_index, :), 'o', ...
            'Color', colors_i_f(i_f_index, :), 'HandleVisibility', 'off');
        
        legend('off');
        legend('show');
    end
    
    title(['Curvas de pot{\^{e}}ncia indexadas pela tens{\~{a}}o da carga ($n_{r}$ $=$ $' num2str(n_r_list(n_r_index)) '$ $[rpm]$)']);
    xlabel('$v_o$ $[V]$');
    ylabel('$P_o$ $[W]$');
    leg = legend;
    leg.Location = 'NorthWest';
    grid on;
end

%% Armazenamento dos resultados de simulação

save('results/PowerAnalysis/P_v_o.mat', 'simIn', 'simOut', 'P_v_o_ana', 'P_v_o_sim', ...
    'P_v_o_mpp_ana', 'P_v_o_mpp_sim', 'v_o_mpp_ana', 'v_o_mpp_sim', '-v7.3');

%% Configuração dos casos de teste para carga de impedância constante

% Configuração dos casos de teste como entrada do modelo analítico
[i_f_grid, n_r_grid, z_o_grid] = meshgrid(i_f_list, n_r_list, z_o_list);

% Configuração dos casos de teste como entrada do modelo no Simulink
load_model = 2;
clear simIn;

for i_f_index = 1:length(i_f_list)
    i_f = i_f_grid(1, i_f_index, 1);
    
    for n_r_index = 1:length(n_r_list)
        n_r = n_r_grid(n_r_index, 1, 1);
        
        for z_o_index = 1:length(z_o_list)
            z_o = z_o_grid(1, 1, z_o_index);
            
            simIn(n_r_index, i_f_index, z_o_index) = Simulink.SimulationInput('PowerAnalysis');
            simIn(n_r_index, i_f_index, z_o_index) = simIn(n_r_index, i_f_index, z_o_index).setVariable('load_model', load_model);
            simIn(n_r_index, i_f_index, z_o_index) = simIn(n_r_index, i_f_index, z_o_index).setBlockParameter('PowerAnalysis/i_f', 'Value', num2str(i_f));
            simIn(n_r_index, i_f_index, z_o_index) = simIn(n_r_index, i_f_index, z_o_index).setBlockParameter('PowerAnalysis/n_r', 'Value', num2str(n_r));
            simIn(n_r_index, i_f_index, z_o_index) = simIn(n_r_index, i_f_index, z_o_index).setBlockParameter('PowerAnalysis/Load/Impedance/z_o', 'R', num2str(z_o));
        end
    end
end

% Transformação de matriz de entradas em vetor
simIn = reshape(simIn, [length(n_r_list)*length(i_f_list)*length(z_o_list) 1]);

%% Análise do efeito da variação da impedância da carga

% Modelo analítico
P_z_o_ana = P_z_o(n_r_grid, i_f_grid, z_o_grid);
P_z_o_ana(imag(P_z_o_ana) ~= 0) = 0;

% Execução da simulação paralelizada
simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
    'TransferBaseWorkspaceVariables', 'on');

simIn = reshape(simIn, [length(n_r_list) length(i_f_list) length(z_o_list)]);
simOut = reshape(simOut, [length(n_r_list) length(i_f_list) length(z_o_list)]);

for i_f_index = 1:length(i_f_list)
    for n_r_index = 1:length(n_r_list)
        for z_o_index = 1:length(z_o_list)
            P_z_o_sim(n_r_index, i_f_index, z_o_index) = mean(simOut(n_r_index, i_f_index, z_o_index).p_l.data(round(end/2):end));
        end
    end
end

%% Identificação dos pontos de máxima potência indexados pela impedância da carga

% 
[P_z_o_mpp_ana, z_o_mpp_index_ana] = max(P_z_o_ana, [], 3);
z_o_mpp_ana = z_o_list(z_o_mpp_index_ana);

% 
[P_z_o_mpp_sim, z_o_mpp_index_sim] = max(P_z_o_sim, [], 3);
z_o_mpp_sim = z_o_list(z_o_mpp_index_sim);

%% Traço dos resultados relativos à variação da impedância da carga

% 
for i_f_index = 1:length(i_f_list)
    
    figure_index = figure_index + 1;
    figure(figure_index)
    
    for n_r_index = 1:length(n_r_list)
        
        plot(z_o_list, squeeze(P_z_o_ana(n_r_index, i_f_index, :)), '-', 'Color', colors_n_r(n_r_index, :), ...
            'DisplayName', ['$n_{r} = ' num2str(n_r_list(n_r_index)) ' rpm$']);
        hold on;
        plot(z_o_mpp_ana(n_r_index, i_f_index, :), P_z_o_mpp_ana(n_r_index, i_f_index, :), 'o', ...
            'Color', colors_n_r(n_r_index, :), 'HandleVisibility', 'off');
        hold on;
        plot(z_o_list, squeeze(P_z_o_sim(n_r_index, i_f_index, :)), '--', 'Color', colors_n_r(n_r_index, :), ...
            'HandleVisibility', 'off');
        hold on;
        plot(z_o_mpp_sim(n_r_index, i_f_index, :), P_z_o_mpp_sim(n_r_index, i_f_index, :), 'o', ...
            'Color', colors_n_r(n_r_index, :), 'HandleVisibility', 'off');
        
        legend('off');
        legend('show');
    end
    
    title(['Curvas de pot{\^{e}}ncia indexadas pela imped{\^{a}}ncia da carga ($i_{f}$ $=$ $' num2str(i_f_list(i_f_index)) '$ $[A]$)']);
    xlabel('$z_o$ $[\Omega]$');
    ylabel('$P_o$ $[W]$');
    leg = legend;
    leg.Location = 'NorthEast';
    grid on;
end

% 
for n_r_index = 1:length(n_r_list)
    
    figure_index = figure_index + 1;
    figure(figure_index)
    
    for i_f_index = 1:length(i_f_list)
        
        plot(z_o_list, squeeze(P_z_o_ana(n_r_index, i_f_index, :)), '-', 'Color', colors_i_f(i_f_index, :), ...
            'DisplayName', ['$i_{f} = ' num2str(i_f_list(i_f_index)) ' A$']);
        hold on;
        plot(z_o_mpp_ana(n_r_index, i_f_index, :), P_z_o_mpp_ana(n_r_index, i_f_index, :), 'o', ...
            'Color', colors_i_f(i_f_index, :), 'HandleVisibility', 'off');
        hold on;
        plot(z_o_list, squeeze(P_z_o_sim(n_r_index, i_f_index, :)), '--', 'Color', colors_i_f(i_f_index, :), ...
            'HandleVisibility', 'off');
        hold on;
        plot(z_o_mpp_sim(n_r_index, i_f_index, :), P_z_o_mpp_sim(n_r_index, i_f_index, :), 'o', ...
            'Color', colors_i_f(i_f_index, :), 'HandleVisibility', 'off');
        
        legend('off');
        legend('show');
    end
    
    title(['Curvas de pot{\^{e}}ncia indexadas pela imped{\^{a}}ncia da carga ($n_{r}$ $=$ $' num2str(n_r_list(n_r_index)) '$ $[rpm]$)']);
    xlabel('$z_o$ $[\Omega]$');
    ylabel('$P_o$ $[W]$');
    leg = legend;
    leg.Location = 'NorthEast';
    grid on;
end

%% Armazenamento dos resultados de simulação

save('results/PowerAnalysis/P_z_o.mat', 'simIn', 'simOut', 'P_z_o_ana', 'P_z_o_sim', ...
    'P_z_o_mpp_ana', 'P_z_o_mpp_sim', 'z_o_mpp_ana', 'z_o_mpp_sim', '-v7.3');

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = sprintf('results/PowerAnalysis/Figura_%d', i);
    saveFigure(figure(i), fileName, 'fig');
end

%% Armazenamento dos resultados de simulação

save('results/PowerAnalysis/simEnv.mat', 'alternator', 'rectifier', ...
    'i_f_list', 'n_r_list', 'v_o_list', 'z_o_list', '-v7.3');
