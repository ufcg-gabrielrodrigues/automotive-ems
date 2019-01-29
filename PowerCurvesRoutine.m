%% Parâmetros temporais

T_s = 1e-6; % Passo de cálculo utilizado pelo 'solver' [s]

%% Varredura de parâmetros

% Lista de parâmetros a serem varridos individualmente
i_f_max = 4.5;                    	% Corrente de excitação máxima [A]
n_r_list = (2000:500:7500)';       	% Velocidade do alternador [rpm]
v_o_list = (0.0:1.0:80.0)';         % Tensão de saída [V]
z_o_list = [0.01 0.05:0.05:2.0]';   % Impedância de saída [Ohm]

%% Alternador

% Fator de acoplamento
if (isfield(alternator.k_e, 'function'))
    k_e = @(i_f) alternator.k_e.function(i_f);
else
    k_e = @(i_f) alternator.k_e.value;
end

if (alternator.stator.connection == delta)
    k_e = @(i_f) k_e(i_f)./sqrt(3);
end

% Indutância de estator
if (isfield(alternator.stator.l, 'function'))
    l_s = @(i_f) alternator.stator.l.function(i_f);
else
    l_s = @(i_f) alternator.stator.l.value;
end

if (alternator.stator.connection == delta)
    l_s = @(i_f) l_s(i_f)./3;
end

%% 

omega_e = @(n_r) n_r.*(2.*pi./60).*alternator.p;
v_s = @(n_r, i_f) k_e(i_f).*omega_e(n_r).*i_f;

%% 

P_v_o = @(n_r, i_f, v_o) (3.*v_o./pi).*(sqrt(v_s(n_r, i_f).^2 - (2.*v_o./pi).^2))./(omega_e(n_r).*l_s(i_f));
P_z_o = @(n_r, i_f, z_o) ((3.*pi.*v_s(n_r, i_f)).^2.*z_o)./((pi.^2.*omega_e(n_r).*l_s(i_f)).^2 + (6.*z_o).^2);

%% 

% 
[n_r, v_o] = meshgrid(n_r_list, v_o_list);

P_v_o = P_v_o(n_r, i_f_max, v_o);
P_v_o(imag(P_v_o) ~= 0) = 0;

% 
figure(1)

for curve_index = 1:length(n_r_list)
    plot(v_o_list, P_v_o(:, curve_index));
    hold on;
end

grid on;

%% 

[n_r, z_o] = meshgrid(n_r_list, z_o_list);

P_z_o = P_z_o(n_r, i_f_max, z_o);
P_z_o(imag(P_z_o) ~= 0) = 0;

% 
figure(2)

for curve_index = 1:length(n_r_list)
    plot(z_o_list, P_z_o(:, curve_index));
    hold on;
end

grid on;

%% Curvas de potência 

% Comparadas por corrente de excitação
for n_r = n_r_list'
    
    for r_l = r_l_list'
        
        colors = distinguishable_colors(length(i_f_list));
        color_index = 0;
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for i_f = i_f_list'
            color_index = color_index + 1;
            curve_indexes = find((test_case_matrix(:, 1) == i_f) ...
                & (test_case_matrix(:, 2) == n_r) ...
                & (test_case_matrix(:, 3) == r_l));
            
            plot(test_case_matrix(curve_indexes, 4), test_case_matrix(curve_indexes, 5), ...
                'Color', colors(color_index, :), 'DisplayName', ...
                ['$i_{f} = ' num2str(i_f) ' A$']);
            hold on;
            
            legend('off');
            legend('show');
        end
        
        title(['Curvas de pot{\^{e}}ncia [$n_{r} = ' num2str(n_r) ...
            ' rpm$; $r_{l} = ' num2str(r_l) ' \Omega$]']);
        xlabel('$u$');
        ylabel('$P [W]$');
        leg = legend;
        leg.Location = 'SouthOutside';
        leg.NumColumns = ceil(length(i_f_list)/leg_entries_per_columns);
        grid on;
    end
end

% Comparadas por velocidade do alternador
for i_f = i_f_list'
    
    for r_l = r_l_list'
        
        colors = distinguishable_colors(length(n_r_list));
        color_index = 0;
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for n_r = n_r_list'
            color_index = color_index + 1;
            curve_indexes = find((test_case_matrix(:, 1) == i_f) ...
                & (test_case_matrix(:, 2) == n_r) ...
                & (test_case_matrix(:, 3) == r_l));
            
            plot(test_case_matrix(curve_indexes, 4), test_case_matrix(curve_indexes, 5), ...
                'Color', colors(color_index, :), 'DisplayName', ...
                ['$n_{r} = ' num2str(n_r) ' rpm$']);
            hold on;
            
            legend('off');
            legend('show');
        end
        
        title(['Curvas de pot{\^{e}}ncia [$i_{f} = ' num2str(i_f) ...
            ' A$; $r_{l} = ' num2str(r_l) ' \Omega$]']);
        xlabel('$u$');
        ylabel('$P [W]$');
        leg = legend;
        leg.Location = 'SouthOutside';
        leg.NumColumns = ceil(length(n_r_list)/leg_entries_per_columns);
        grid on;
    end
end

% Comparadas por carga
for i_f = i_f_list'
    
    for n_r = n_r_list'
        
        colors = distinguishable_colors(length(r_l_list));
        color_index = 0;
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        for r_l = r_l_list'
            color_index = color_index + 1;
            curve_indexes = find((test_case_matrix(:, 1) == i_f) ...
                & (test_case_matrix(:, 2) == n_r) ...
                & (test_case_matrix(:, 3) == r_l));
            
            plot(test_case_matrix(curve_indexes, 4), test_case_matrix(curve_indexes, 5), ...
                'Color', colors(color_index, :), 'DisplayName', ...
                ['$r_{l} = ' num2str(r_l, '%1.2f') ' \Omega$']);
            hold on;
            
            legend('off');
            legend('show');
        end
        
        title(['Curvas de pot{\^{e}}ncia [$i_{f} = ' num2str(i_f) ...
            ' A$; $n_{r} = ' num2str(n_r) ' rpm$]']);
        xlabel('$u$');
        ylabel('$P [W]$');
        leg = legend;
        leg.Location = 'SouthOutside';
        leg.NumColumns = ceil(length(r_l_list)/leg_entries_per_columns);
        grid on;
    end
end

%% Relação dos pontos de máxima potência com variáveis do sistema

% Relação com a corrente de excitação
for n_r = n_r_list'
    
    for r_l = r_l_list'
        
        curve_indexes = find((mpp_matrix(:, 2) == n_r) & (mpp_matrix(:, 3) == r_l));
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        subplot(2, 1, 1)
        plot(mpp_matrix(curve_indexes, 1), mpp_matrix(curve_indexes, 4), 'o-');
        xlabel('$i_{f} [A]$');
        ylabel('$u$');
        ylim([0 1]);
        grid on;
        
        subplot(2, 1, 2)
        plot(mpp_matrix(curve_indexes, 1), mpp_matrix(curve_indexes, 5), 'o-');
        xlabel('$i_{f} [A]$');
        ylabel('$P [W]$');
        grid on;
        
        suptitle(['Rela{\c{c}}{\~{a}}o entre o ponto de m{\''{a}}xima ' ...
            'pot{\^{e}}ncia e a corrente de excita{\c{c}}{\~{a}}o [$n_{r} = ' ...
            num2str(n_r) ' rpm$; $r_{l} = ' num2str(r_l) ' \Omega$]']);
    end
end

% Relação com a velocidade do alternador
for i_f = i_f_list'
    
    for r_l = r_l_list'
        
        curve_indexes = find((mpp_matrix(:, 1) == i_f) & (mpp_matrix(:, 3) == r_l));
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        subplot(2, 1, 1)
        plot(mpp_matrix(curve_indexes, 2), mpp_matrix(curve_indexes, 4), 'o-');
        xlabel('$n_{r} [rpm]$');
        ylabel('$u$');
        ylim([0 1]);
        grid on;
        
        subplot(2, 1, 2)
        plot(mpp_matrix(curve_indexes, 2), mpp_matrix(curve_indexes, 5), 'o-');
        xlabel('$n_{r} [rpm]$');
        ylabel('$P [W]$');
        grid on;
        
        suptitle(['Rela{\c{c}}{\~{a}}o entre o ponto de m{\''{a}}xima ' ...
            'pot{\^{e}}ncia e a velocidade do alternador [$i_{f} = ' ...
            num2str(i_f) ' A$; $r_{l} = ' num2str(r_l) ' \Omega$]']);
    end
end

% Relação com a impedância de carga
for i_f = i_f_list'
    
    for n_r = n_r_list'
        
        curve_indexes = find((mpp_matrix(:, 1) == i_f) & (mpp_matrix(:, 2) == n_r));
        
        figure_index = figure_index + 1;
        figure(figure_index)
        
        subplot(2, 1, 1)
        plot(mpp_matrix(curve_indexes, 3), mpp_matrix(curve_indexes, 4), 'o-');
        xlabel('$r_{l} [\Omega]$');
        ylabel('$u$');
        ylim([0 1]);
        grid on;
        
        subplot(2, 1, 2)
        plot(mpp_matrix(curve_indexes, 3), mpp_matrix(curve_indexes, 5), 'o-');
        xlabel('$r_{l} [\Omega]$');
        ylabel('$P [W]$');
        grid on;
        
        suptitle(['Rela{\c{c}}{\~{a}}o entre o ponto de m{\''{a}}xima ' ...
            'pot{\^{e}}ncia e a imped{\^{a}}ncia de carga [$i_{f} = ' ...
            num2str(i_f) ' A$; $n_{r} = ' num2str(n_r) ' rpm$]']);
    end
end

%% Armazenamento de figuras

for i = 1:figure_index
    fileName = sprintf('results/MPPTCurves/Figura_%d', i);
    saveFigure(figure(i), fileName, 'fig');
end

%% Armazenamento dos resultados de simulação

% Tratamento para simulação particionada
if (~sim_split_flag)
    save('results/MPPTCurves/test_case_out.mat', 'test_case_out', '-v7.3');
end

save('results/MPPTCurves/test_case_matrix.mat', 'test_case_matrix', '-v7.3');
save('results/MPPTCurves/mpp_matrix.mat', 'mpp_matrix', '-v7.3');
save('results/MPPTCurves/mpp_map.mat', 'i_f_list', 'n_r_list', 'r_l_list', ...
    'mpp_u_3d', 'mpp_p_3d', '-v7.3');
