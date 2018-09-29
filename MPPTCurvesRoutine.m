%%

i_f_list = [1; 3];                      % [A]
n_alt_list = [1800; 3000; 4000; 5000];  % [rpm]
u_list = linspace(0, 1, 21);            % []

% i_f_list = [1];           	% [A]
% n_alt_list = [1800];        % [rpm]
% u_list = linspace(0, 1, 2);	% []

%%

param_sweep = [];

for index_i_f = 1:length(i_f_list)
    
    n_alt_sweep = [];
    
    for index_n_alt = 1:length(n_alt_list)
        n_alt_tmp = n_alt_list(index_n_alt) * ones(length(u_list), 1);
        n_alt_tmp = [n_alt_tmp, u_list'];
        n_alt_sweep = [n_alt_sweep; n_alt_tmp];
    end
    
    i_f_tmp = i_f_list(index_i_f) * ones(length(n_alt_sweep), 1);
    i_f_tmp = [i_f_tmp, n_alt_sweep];
    param_sweep = [param_sweep; i_f_tmp];
end

%% Inicializa modelo no Simulink

open_system('MPPTCurves.slx', 'loadonly');

%% Parâmetros de simulação

% Parâmetros do 'solver' local para sistemas físicos
T_s = 1e-6;
set_param('MPPTCurves/Solver Configuration', 'UseLocalSolver', 'on');
set_param('MPPTCurves/Solver Configuration', 'LocalSolverChoice', 'NE_TRAPEZOIDAL_ADVANCER');
set_param('MPPTCurves/Solver Configuration', 'LocalSolverSampleTime', num2str(T_s));
set_param('MPPTCurves/Solver Configuration', 'DoFixedCost', 'on');
set_param('MPPTCurves/Solver Configuration', 'MaxNonlinIter', '20');

% Parâmetros do 'solver' global
simulationParameters.StopTime = '1e-1'; % [s]

%% Execução da simulação em ambiente Simulink

[num_cases, ~] = size(param_sweep);

for test_case_index = 1:num_cases
    test_case = param_sweep(test_case_index, :);
    i_f = test_case(1);
    n_alt = test_case(2);
    u = test_case(3);
    
    %
    simout{test_case_index} = sim('MPPTCurves', simulationParameters);
end

%% Salva e finaliza modelo no Simulink

save_system('MPPTCurves.slx');
close_system('MPPTCurves.slx');

%% Registro de resultados obtidos no caso de teste

for test_case_index = 1:num_cases
    % 
    test_case = param_sweep(test_case_index, :);
    i_f = test_case(1);
    n_alt = test_case(2);
    u = test_case(3);
    
    % Alternador
    test_case_out(test_case_index).alternator.rotor.n = n_alt;
    test_case_out(test_case_index).alternator.rotor.l.i = i_f;
    
    test_case_out(test_case_index).alternator.stator.input.e.value = simout{test_case_index}.e_a_abc;
    test_case_out(test_case_index).alternator.stator.output.v = simout{test_case_index}.v_a_abc;
    test_case_out(test_case_index).alternator.stator.output.i = simout{test_case_index}.i_a_abc;
    
    % Retificador
    test_case_out(test_case_index).rectifier.control.u = u;
    
    % Carga
    test_case_out(test_case_index).load.v.curve = simout{test_case_index}.v_l;
    test_case_out(test_case_index).load.i.curve = simout{test_case_index}.i_l;
    test_case_out(test_case_index).load.p.curve = simout{test_case_index}.p_l;
    
    num_samples = length(test_case_out(test_case_index).load.p.curve.Data);
    test_case_out(test_case_index).load.v.avg = mean(test_case_out(test_case_index).load.v.curve.Data(floor(num_samples/2 + 1):num_samples));
    test_case_out(test_case_index).load.i.avg = mean(test_case_out(test_case_index).load.i.curve.Data(floor(num_samples/2 + 1):num_samples));
    test_case_out(test_case_index).load.p.avg = mean(test_case_out(test_case_index).load.p.curve.Data(floor(num_samples/2 + 1):num_samples));
end

%% 

i_f_lim = [];
n_alt_lim = [];

for i_f_index = 1:length(i_f_list)
    % 
    i_f_lim_case = find(param_sweep(:, 1) == i_f_list(i_f_index));
    i_f_lim(i_f_index, :) = [i_f_lim_case(1) i_f_lim_case(end)];
    
    % 
    for n_alt_index = 1:length(n_alt_list)
        %
        n_alt_lim_case = find(param_sweep(i_f_lim(i_f_index, 1):i_f_lim(i_f_index, 2), 2) == n_alt_list(n_alt_index));
        n_alt_lim(n_alt_index + (i_f_index - 1)*length(n_alt_list), :) = [n_alt_lim_case(1) n_alt_lim_case(end)] + (i_f_lim(i_f_index, 1) - 1);
    end
end

curve_lim = n_alt_lim;

%% 

for curve_index = 1:length(curve_lim)
    mppt_curve(curve_index).n_alt = test_case_out(curve_lim(curve_index, 1)).alternator.rotor.n;
    mppt_curve(curve_index).i_f = test_case_out(curve_lim(curve_index, 1)).alternator.rotor.l.i;
    
    for point_index = curve_lim(curve_index, 1):curve_lim(curve_index, 2)
        mppt_curve(curve_index).trace(point_index - curve_lim(curve_index, 1) + 1, 1) = test_case_out(point_index).rectifier.control.u;
        mppt_curve(curve_index).trace(point_index - curve_lim(curve_index, 1) + 1, 2) = test_case_out(point_index).load.p.avg;
    end
    
    % 
    figure(curve_index)
    plot(mppt_curve(curve_index).trace(:,1), mppt_curve(curve_index).trace(:,2));
    title(['Curva P x u (i_{f} = ' num2str(mppt_curve(curve_index).i_f) ...
        ' A; n_{alt} = ' num2str(mppt_curve(curve_index).n_alt) ' rpm)']);
    xlabel('u');
    ylabel('P_{l}');
    grid on;
end

%% Armazenamento dos resultados de simulação

save('results/test_case_out.mat', 'test_case_out', '-v7.3');
save('results/mppt_curve.mat', 'mppt_curve', '-v7.3');
