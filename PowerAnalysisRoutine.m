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
alternator.stator.r.T = 150;    % [oC]

% Indutância mútua
if (isfield(alternator.m_f, 'function'))
    m_f_fun = @(i_f) alternator.m_f.function(i_f);
else
    m_f_fun = @(i_f) alternator.m_f.value;
end

if (alternator.stator.connection == delta)
    m_f_fun = @(i_f) m_f_fun(i_f)./sqrt(3);
end

% Indutância por fase da armadura
if (isfield(alternator.stator.l, 'function'))
    l_s_fun = @(i_f) alternator.stator.l.function(i_f);
else
    l_s_fun = @(i_f) alternator.stator.l.value;
end

if (alternator.stator.connection == delta)
    l_s_fun = @(i_f) l_s_fun(i_f)./3;
end

% Função de cálculo da frequência elétrica
omega = @(n_r) n_r.*(2.*pi./60).*alternator.p;

% Função de cálculo da tensão induzida no estator
v_s = @(n_r, i_f) m_f_fun(i_f).*omega(n_r).*i_f;

%% Modelos analíticos para cálculo de potência

% Carga de tensão constante
p_v_dc = @(n_r, i_f, v_dc) (3.*v_dc./pi).*(sqrt(v_s(n_r, i_f).^2 - (2.*v_dc./pi).^2))./(omega(n_r).*l_s_fun(i_f));

% Carga de impedância constante
p_z_dc = @(n_r, i_f, z_dc) ((3.*pi.*v_s(n_r, i_f)).^2.*z_dc)./((pi.^2.*omega(n_r).*l_s_fun(i_f)).^2 + (6.*z_dc).^2);

%% Varredura de parâmetros

% Lista de parâmetros a serem varridos individualmente
i_f_list = [0.01 1.0:1.0:5.0]';     % Corrente de excitação máxima [A]
n_r_list = (2000:500:7500)';        % Velocidade do alternador [rpm]
v_dc_list = (0.0:1.0:80.0)';        % Tensão de saída [V]
z_dc_list = [0.01 0.05:0.05:2.0]';	% Impedância de saída [Ohm]

%% Parâmetros auxiliares para figuras

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

%% Análise para carga tipo tensão

if (~isempty(v_dc_list))
    %% Configuração dos casos de teste para carga de tensão constante
    
    % Configuração dos casos de teste como entrada do modelo analítico
    [i_f_grid, n_r_grid, v_dc_grid] = meshgrid(i_f_list, n_r_list, v_dc_list);
    
    % Configuração dos casos de teste como entrada do modelo no Simulink
    load_model = 1;
    clear simIn;
    
    for i_f_index = 1:length(i_f_list)
        i_f = i_f_grid(1, i_f_index, 1);
        
        for n_r_index = 1:length(n_r_list)
            n_r = n_r_grid(n_r_index, 1, 1);
            
            for v_dc_index = 1:length(v_dc_list)
                v_dc = v_dc_grid(1, 1, v_dc_index);
                
                simIn(n_r_index, i_f_index, v_dc_index) = Simulink.SimulationInput('PowerAnalysis');
                simIn(n_r_index, i_f_index, v_dc_index) = simIn(n_r_index, i_f_index, v_dc_index).setVariable('load_model', load_model);
                simIn(n_r_index, i_f_index, v_dc_index) = simIn(n_r_index, i_f_index, v_dc_index).setBlockParameter('PowerAnalysis/i_f', 'Value', num2str(i_f));
                simIn(n_r_index, i_f_index, v_dc_index) = simIn(n_r_index, i_f_index, v_dc_index).setBlockParameter('PowerAnalysis/n_r', 'Value', num2str(n_r));
                simIn(n_r_index, i_f_index, v_dc_index) = simIn(n_r_index, i_f_index, v_dc_index).setBlockParameter('PowerAnalysis/Load/Voltage/v_dc', 'DC', num2str(v_dc));
            end
        end
    end
    
    % Transformação de matriz de entradas em vetor
    simIn = reshape(simIn, [length(n_r_list)*length(i_f_list)*length(v_dc_list) 1]);
    
    %% Análise do efeito da variação da tensão na carga
    
    % Modelo analítico
    p_v_dc_ana = p_v_dc(n_r_grid, i_f_grid, v_dc_grid);
    p_v_dc_ana(imag(p_v_dc_ana) ~= 0) = 0;
    
    % Execução da simulação paralelizada
    simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
        'TransferBaseWorkspaceVariables', 'on');
    
    simIn = reshape(simIn, [length(n_r_list) length(i_f_list) length(v_dc_list)]);
    simOut = reshape(simOut, [length(n_r_list) length(i_f_list) length(v_dc_list)]);
    
    for i_f_index = 1:length(i_f_list)
        for n_r_index = 1:length(n_r_list)
            for v_dc_index = 1:length(v_dc_list)
                if (isempty(simOut(n_r_index, i_f_index, v_dc_index).ErrorMessage))
                    p_v_dc_sim(n_r_index, i_f_index, v_dc_index) = mean(simOut(n_r_index, i_f_index, v_dc_index).p_dc.data(round(end/2):end));
                else
                    p_v_dc_sim(n_r_index, i_f_index, v_dc_index) = nan;
                end
            end
        end
    end
    
    %% Identificação dos pontos de máxima potência indexados pela tensão da carga
    
    %
    [p_v_dc_mpp_ana, v_dc_mpp_index_ana] = max(p_v_dc_ana, [], 3);
    v_dc_mpp_ana = v_dc_list(v_dc_mpp_index_ana);
    
    %
    [p_v_dc_mpp_sim, v_dc_mpp_index_sim] = max(p_v_dc_sim, [], 3);
    v_dc_mpp_sim = v_dc_list(v_dc_mpp_index_sim);
    
    %% Traço dos resultados relativos à variação da tensão da carga
    
    % �?ndice de figuras
    figure_index = 0;
    
    % Lista de t�tulos
    tit = cell(length(i_f_list) + length(n_r_list) + 2);
    
    %
    for i_f_index = 1:length(i_f_list)
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for n_r_index = 1:length(n_r_list)
            
            plot(v_dc_list, squeeze(p_v_dc_ana(n_r_index, i_f_index, :)), '-', 'Color', colors_n_r(n_r_index, :), ...
                'DisplayName', ['$n_{r} = ' num2str(n_r_list(n_r_index)) '\,\textrm{rpm}$']);
            hold on;
            plot(v_dc_mpp_ana(n_r_index, i_f_index, :), p_v_dc_mpp_ana(n_r_index, i_f_index, :), 'o', ...
                'Color', colors_n_r(n_r_index, :), 'HandleVisibility', 'off');
            hold on;
            plot(v_dc_list, squeeze(p_v_dc_sim(n_r_index, i_f_index, :)), '--', 'Color', colors_n_r(n_r_index, :), ...
                'HandleVisibility', 'off');
            hold on;
            plot(v_dc_mpp_sim(n_r_index, i_f_index, :), p_v_dc_mpp_sim(n_r_index, i_f_index, :), 'o', ...
                'Color', colors_n_r(n_r_index, :), 'HandleVisibility', 'off');
            
            legend('off');
            legend('show');
        end
        
        tit{figure_index} = ['mpp-vo-if-' num2str(i_f_list(i_f_index))];
        xlabel('$v_{dc}\,[\textrm{V}]$');
        ylabel('$p_{dc}\,[\textrm{W}]$');
        leg = legend;
        leg.Location = 'NorthWest';
        grid on;
    end
    
    %
    for n_r_index = 1:length(n_r_list)
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for i_f_index = 1:length(i_f_list)
            
            plot(v_dc_list, squeeze(p_v_dc_ana(n_r_index, i_f_index, :)), '-', 'Color', colors_i_f(i_f_index, :), ...
                'DisplayName', ['$i_{f} = ' num2str(i_f_list(i_f_index), '%1.2f') '\,\textrm{A}$']);
            hold on;
            plot(v_dc_mpp_ana(n_r_index, i_f_index, :), p_v_dc_mpp_ana(n_r_index, i_f_index, :), 'o', ...
                'Color', colors_i_f(i_f_index, :), 'HandleVisibility', 'off');
            hold on;
            plot(v_dc_list, squeeze(p_v_dc_sim(n_r_index, i_f_index, :)), '--', 'Color', colors_i_f(i_f_index, :), ...
                'HandleVisibility', 'off');
            hold on;
            plot(v_dc_mpp_sim(n_r_index, i_f_index, :), p_v_dc_mpp_sim(n_r_index, i_f_index, :), 'o', ...
                'Color', colors_i_f(i_f_index, :), 'HandleVisibility', 'off');
            
            legend('off');
            legend('show');
        end
        
        tit{figure_index} = ['mpp-vo-nr-' num2str(n_r_list(n_r_index))];
        xlabel('$v_{dc}\,[\textrm{V}]$');
        ylabel('$p_{dc}\,[\textrm{W}]$');
        leg = legend;
        leg.Location = 'NorthWest';
        grid on;
    end
    
    %
    try
        figure_index = figure_index + 1;
        figure(figure_index)
        
        h_mpp_ana = surf(squeeze(n_r_grid(:, :, 1)), squeeze(i_f_grid(:, :, 1)), v_dc_mpp_ana);
        colormap(spring);
        freezeColors;
        hold on;
        h_mpp_sim = surf(squeeze(n_r_grid(:, :, 1)), squeeze(i_f_grid(:, :, 1)), v_dc_mpp_sim);
        colormap(winter);
        
        tit{figure_index} = 'mpp-vo-sup';
        xlabel('$n_{r}\,[\textrm{rpm}]$');
        ylabel('$i_{f}\,[\textrm{A}]$');
        zlabel('$v_{dc}\,[\textrm{V}]$');
        legend([h_mpp_ana, h_mpp_sim], {'Superf{\''i}cie obtida analiticamente', ...
            'Superf{\''i}cie obtida via simula{\c{c}}{\~{a}}o'}, 'Location', 'NorthEast');
        grid on;
    catch
        close;
        figure_index = figure_index - 1;
    end
    
    try
        figure_index = figure_index + 1;
        figure(figure_index)
        
        [xData, yData, zData] = prepareSurfaceData(n_r_list, i_f_list, p_v_dc_mpp_ana);
        [p_v_dc_mpp_ana_fit, ~] = fit([xData, yData], zData, 'thinplateinterp', 'Normalize', 'on');
        
        [xData, yData, zData] = prepareSurfaceData(n_r_list, i_f_list, p_v_dc_mpp_sim);
        [p_v_dc_mpp_sim_fit, ~] = fit([xData, yData], zData, 'thinplateinterp', 'Normalize', 'on');
        
        h_mpp_ana = surf(squeeze(n_r_grid(:, :, 1)), squeeze(i_f_grid(:, :, 1)), p_v_dc_mpp_ana);
        colormap(spring);
        
        freezeColors;
        hold on;
        
        h_mpp_sim = surf(squeeze(n_r_grid(:, :, 1)), squeeze(i_f_grid(:, :, 1)), p_v_dc_mpp_sim);
        colormap(winter);
        
        tit{figure_index} = 'mpp-po-v-sup';
        xlabel('$n_{r}\,[\textrm{rpm}]$');
        ylabel('$i_{f}\,[\textrm{A}]$');
        zlabel('$p_{dc}\,[\textrm{W}]$');
        legend([h_mpp_ana, h_mpp_sim], {'Superf{\''i}cie obtida analiticamente', ...
            'Superf{\''i}cie obtida via simula{\c{c}}{\~{a}}o'}, 'Location', 'NorthEast');
        grid on;
    catch
        close
        figure_index = figure_index - 1;
    end
    
    %% Armazenamento de figuras
    
    for i = 1:figure_index
        fileName = ['results/PowerAnalysis/' tit{i}];
        saveFigure(figure(i), fileName, 'fig');
        saveFigure(figure(i), fileName, 'eps');
    end
    
    close all;
    
    %% Armazenamento dos resultados de simulação
    
    save('results/PowerAnalysis/p_v_dc.mat', 'simIn', 'simOut', 'p_v_dc_ana', 'p_v_dc_sim', ...
        'p_v_dc_mpp_ana', 'p_v_dc_mpp_sim', 'v_dc_mpp_ana', 'v_dc_mpp_sim', '-v7.3');
end

%% Análise para carga tipo impedância

if (~isempty(z_dc_list))
    %% Configuração dos casos de teste para carga de impedância constante
    
    % Configuração dos casos de teste como entrada do modelo analítico
    [i_f_grid, n_r_grid, z_dc_grid] = meshgrid(i_f_list, n_r_list, z_dc_list);
    
    % Configuração dos casos de teste como entrada do modelo no Simulink
    load_model = 2;
    clear simIn;
    
    for i_f_index = 1:length(i_f_list)
        i_f = i_f_grid(1, i_f_index, 1);
        
        for n_r_index = 1:length(n_r_list)
            n_r = n_r_grid(n_r_index, 1, 1);
            
            for z_dc_index = 1:length(z_dc_list)
                z_dc = z_dc_grid(1, 1, z_dc_index);
                
                simIn(n_r_index, i_f_index, z_dc_index) = Simulink.SimulationInput('PowerAnalysis');
                simIn(n_r_index, i_f_index, z_dc_index) = simIn(n_r_index, i_f_index, z_dc_index).setVariable('load_model', load_model);
                simIn(n_r_index, i_f_index, z_dc_index) = simIn(n_r_index, i_f_index, z_dc_index).setBlockParameter('PowerAnalysis/i_f', 'Value', num2str(i_f));
                simIn(n_r_index, i_f_index, z_dc_index) = simIn(n_r_index, i_f_index, z_dc_index).setBlockParameter('PowerAnalysis/n_r', 'Value', num2str(n_r));
                simIn(n_r_index, i_f_index, z_dc_index) = simIn(n_r_index, i_f_index, z_dc_index).setBlockParameter('PowerAnalysis/Load/Impedance/z_dc', 'R', num2str(z_dc));
            end
        end
    end
    
    % Transformação de matriz de entradas em vetor
    simIn = reshape(simIn, [length(n_r_list)*length(i_f_list)*length(z_dc_list) 1]);
    
    %% Análise do efeito da variação da impedância da carga
    
    % Modelo analítico
    p_z_dc_ana = p_z_dc(n_r_grid, i_f_grid, z_dc_grid);
    p_z_dc_ana(imag(p_z_dc_ana) ~= 0) = 0;
    
    % Execução da simulação paralelizada
    simOut = parsim(simIn, 'ShowProgress', 'on', 'ShowSimulationManager', 'on', ...
        'TransferBaseWorkspaceVariables', 'on');
    
    simIn = reshape(simIn, [length(n_r_list) length(i_f_list) length(z_dc_list)]);
    simOut = reshape(simOut, [length(n_r_list) length(i_f_list) length(z_dc_list)]);
    
    for i_f_index = 1:length(i_f_list)
        for n_r_index = 1:length(n_r_list)
            for z_dc_index = 1:length(z_dc_list)
                if (isempty(simOut(n_r_index, i_f_index, z_dc_index).ErrorMessage))
                    p_z_dc_sim(n_r_index, i_f_index, z_dc_index) = mean(simOut(n_r_index, i_f_index, z_dc_index).p_dc.data(round(end/2):end));
                else
                    p_z_dc_sim(n_r_index, i_f_index, z_dc_index) = nan;
                end
            end
        end
    end
    
    %% Identificação dos pontos de máxima potência indexados pela impedância da carga
    
    %
    [p_z_dc_mpp_ana, z_dc_mpp_index_ana] = max(p_z_dc_ana, [], 3);
    z_dc_mpp_ana = z_dc_list(z_dc_mpp_index_ana);
    
    %
    [p_z_dc_mpp_sim, z_dc_mpp_index_sim] = max(p_z_dc_sim, [], 3);
    z_dc_mpp_sim = z_dc_list(z_dc_mpp_index_sim);
    
    %% Traço dos resultados relativos à variação da impedância da carga
    
    % �?ndice de figuras
    figure_index = 0;
    
    % Lista de t�tulos
    tit = cell(length(i_f_list) + length(n_r_list) + 2);
    
    %
    for i_f_index = 1:length(i_f_list)
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for n_r_index = 1:length(n_r_list)
            
            plot(z_dc_list, squeeze(p_z_dc_ana(n_r_index, i_f_index, :)), '-', 'Color', colors_n_r(n_r_index, :), ...
                'DisplayName', ['$n_{r} = ' num2str(n_r_list(n_r_index)) '\,\textrm{rpm}$']);
            hold on;
            plot(z_dc_mpp_ana(n_r_index, i_f_index, :), p_z_dc_mpp_ana(n_r_index, i_f_index, :), 'o', ...
                'Color', colors_n_r(n_r_index, :), 'HandleVisibility', 'off');
            hold on;
            plot(z_dc_list, squeeze(p_z_dc_sim(n_r_index, i_f_index, :)), '--', 'Color', colors_n_r(n_r_index, :), ...
                'HandleVisibility', 'off');
            hold on;
            plot(z_dc_mpp_sim(n_r_index, i_f_index, :), p_z_dc_mpp_sim(n_r_index, i_f_index, :), 'o', ...
                'Color', colors_n_r(n_r_index, :), 'HandleVisibility', 'off');
            
            legend('off');
            legend('show');
        end
        
        tit{figure_index} = ['mpp-zo-if-' num2str(i_f_list(i_f_index))];
        xlabel('$z_{dc}\,[\mathrm{\Omega}]$');
        ylabel('$p_{dc}\,[\textrm{W}]$');
        leg = legend;
        leg.Location = 'NorthEast';
        grid on;
    end
    
    %
    for n_r_index = 1:length(n_r_list)
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for i_f_index = 1:length(i_f_list)
            
            plot(z_dc_list, squeeze(p_z_dc_ana(n_r_index, i_f_index, :)), '-', 'Color', colors_i_f(i_f_index, :), ...
                'DisplayName', ['$i_{f} = ' num2str(i_f_list(i_f_index), '%1.2f') '\,\textrm{A}$']);
            hold on;
            plot(z_dc_mpp_ana(n_r_index, i_f_index, :), p_z_dc_mpp_ana(n_r_index, i_f_index, :), 'o', ...
                'Color', colors_i_f(i_f_index, :), 'HandleVisibility', 'off');
            hold on;
            plot(z_dc_list, squeeze(p_z_dc_sim(n_r_index, i_f_index, :)), '--', 'Color', colors_i_f(i_f_index, :), ...
                'HandleVisibility', 'off');
            hold on;
            plot(z_dc_mpp_sim(n_r_index, i_f_index, :), p_z_dc_mpp_sim(n_r_index, i_f_index, :), 'o', ...
                'Color', colors_i_f(i_f_index, :), 'HandleVisibility', 'off');
            
            legend('off');
            legend('show');
        end
        
        tit{figure_index} = ['mpp-zo-nr-' num2str(n_r_list(n_r_index))];
        xlabel('$z_{dc}\,[\mathrm{\Omega}]$');
        ylabel('$p_{dc}\,[\textrm{W}]$');
        leg = legend;
        leg.Location = 'NorthEast';
        grid on;
    end
    
    %
    try
        figure_index = figure_index + 1;
        figure(figure_index)
        
        h_mpp_ana = surf(squeeze(n_r_grid(:, :, 1)), squeeze(i_f_grid(:, :, 1)), z_dc_mpp_ana);
        colormap(spring);
        freezeColors;
        hold on;
        h_mpp_sim = surf(squeeze(n_r_grid(:, :, 1)), squeeze(i_f_grid(:, :, 1)), z_dc_mpp_sim);
        colormap(winter);
        
        tit{figure_index} = 'mpp-zo-sup';
        xlabel('$n_{r}\,[\textrm{rpm}]$');
        ylabel('$i_{f}\,[\textrm{A}]$');
        zlabel('$z_{dc}\,[\mathrm{\Omega}]$');
        legend([h_mpp_ana, h_mpp_sim], {'Superf{\''i}cie obtida analiticamente', ...
            'Superf{\''i}cie obtida via simula{\c{c}}{\~{a}}o'}, 'Location', 'NorthEast');
        grid on;
    catch
        close
        figure_index = figure_index - 1;
    end
    
    try
        figure_index = figure_index + 1;
        figure(figure_index)
        
        [xData, yData, zData] = prepareSurfaceData(n_r_list, i_f_list, p_z_dc_mpp_ana);
        [p_z_dc_mpp_ana_fit, ~] = fit([xData, yData], zData, 'thinplateinterp', 'Normalize', 'on');
        
        [xData, yData, zData] = prepareSurfaceData(n_r_list, i_f_list, p_z_dc_mpp_sim);
        [p_z_dc_mpp_sim_fit, ~] = fit([xData, yData], zData, 'thinplateinterp', 'Normalize', 'on');
        
        h_mpp_ana = surf(squeeze(n_r_grid(:, :, 1)), squeeze(i_f_grid(:, :, 1)), p_z_dc_mpp_ana);
        colormap(spring);
        freezeColors;
        hold on;
        h_mpp_sim = surf(squeeze(n_r_grid(:, :, 1)), squeeze(i_f_grid(:, :, 1)), p_z_dc_mpp_sim);
        colormap(winter);
        
        tit{figure_index} = 'mpp-po-z-sup';
        xlabel('$n_{r}\,[\textrm{rpm}]$');
        ylabel('$i_{f}\,[\textrm{A}]$');
        zlabel('$p_{dc}\,[\textrm{W}]$');
        legend([h_mpp_ana, h_mpp_sim], {'Superf{\''i}cie obtida analiticamente', ...
            'Superf{\''i}cie obtida via simula{\c{c}}{\~{a}}o'}, 'Location', 'NorthEast');
        grid on;
    catch
        close
        figure_index = figure_index - 1;
    end
    
    %% Armazenamento de figuras
    
    for i = 1:figure_index
        fileName = ['results/PowerAnalysis/' tit{i}];
        saveFigure(figure(i), fileName, 'fig');
        saveFigure(figure(i), fileName, 'eps');
    end
    
    close all;
    
    %% Armazenamento dos resultados de simulação
    
    save('results/PowerAnalysis/p_z_dc.mat', 'simIn', 'simOut', 'p_z_dc_ana', 'p_z_dc_sim', ...
        'p_z_dc_mpp_ana', 'p_z_dc_mpp_sim', 'z_dc_mpp_ana', 'z_dc_mpp_sim', 'z_dc_mpp_fit', '-v7.3');
    
end

%% Armazenamento dos resultados de simulação

save('results/PowerAnalysis/simEnv.mat', 'alternator', 'rectifier', ...
    'i_f_list', 'n_r_list', 'v_dc_list', 'z_dc_list', '-v7.3');
