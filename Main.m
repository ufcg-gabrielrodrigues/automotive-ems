%% Script principal da simulação do Sistema de Gerenciamento de Energia (EMS)

%% Preparação de ambiente para simulação

clear all;
close all;
clc;

%% Identificação e adição do diretório principal e seus subdiretórios aos
%  caminhos de busca do MATLAB

% Caminho absoluto do diretório
[root, ~, ~] = fileparts(mfilename('fullpath'));

% Caminho contendo pastas e subpastas
root = genpath(root);

% Inclusão do diretório aos caminhos de busca do MATLAB
addpath(root);

%% Inicializa modelo no Simulink

open_system('AutomotiveEMS.slx', 'loadonly');

%% Motor a combustão interna

% 

%% Alternador

% Determina realização ou não de nova iteração para determinação de
% parâmetros
alternatorFittingFlag = false;

% Escolha do alternador a ser utilizado
alternatorCase = 'Sarafianos2015';

% Executa script que determina parâmetros do alternador
AlternatorParametersEMS;

% Atribui modelo parametrizado da FEM induzida por fase ao bloco
% responsável pelo seu cálculo no ambiente do Simulink
matlabFunctionBlock('AutomotiveEMS/Alternator/Back EMF (abc reference)/backEMF', ...
    sym(alternator.stator.input.e.function), 'FunctionName', 'backEMF', 'outputs', {'e_a'});

% 
matlabFunctionBlock('AutomotiveEMS/Alternator/Stator Inductance (abc reference)/statorInductance', ...
    sym(alternator.stator.l.function), 'FunctionName', 'statorInductance', 'outputs', {'l_s'});

%% Retificador

% 
RectifierParametersEMS;

%% 

iceToAltRotRatio = 2.5;

%% Salva simulação com últimas alterações

save_system('AutomotiveEMS.slx');

%% 

t = 1e+1;           % [s]
solver = 'ode15s';  % 

%% 

simulationParameters.StopTime   = num2str(t);
simulationParameters.Solver     = solver;

%% 

sim('AutomotiveEMS', simulationParameters);

%% 

% Alternador
alternator.rotor.control.q = ans.q_s_f;
alternator.rotor.d.v = timeseries(ans.d_f.Data(:, 1), ans.d_f.Time);
alternator.rotor.d.i = timeseries(ans.d_f.Data(:, 2), ans.d_f.Time);
alternator.rotor.s.v_gs = timeseries(ans.s_f.Data(:, 1), ans.s_f.Time);
alternator.rotor.s.v_ds = timeseries(ans.s_f.Data(:, 2), ans.s_f.Time);
alternator.rotor.s.i = timeseries(ans.s_f.Data(:, 3), ans.s_f.Time);
alternator.rotor.l.i = ans.i_f;

alternator.stator.input.e.value = ans.e_a_abc;
alternator.stator.l.value = ans.l_a_abc;
alternator.stator.l.v = ans.v_l_a_abc;
alternator.stator.l.i = ans.i_a_abc;
alternator.stator.r.v = ans.v_r_a_abc;
alternator.stator.r.i = ans.i_a_abc;
alternator.stator.output.v = ans.v_a_abc;

% Retificador
rectifier.d_1.v = timeseries(ans.d_r_1.Data(:, 1), ans.d_r_1.Time);
rectifier.d_1.i = timeseries(ans.d_r_1.Data(:, 2), ans.d_r_1.Time);
rectifier.d_2.v = timeseries(ans.d_r_2.Data(:, 1), ans.d_r_2.Time);
rectifier.d_2.i = timeseries(ans.d_r_2.Data(:, 2), ans.d_r_2.Time);
rectifier.d_3.v = timeseries(ans.d_r_3.Data(:, 1), ans.d_r_3.Time);
rectifier.d_3.i = timeseries(ans.d_r_3.Data(:, 2), ans.d_r_3.Time);

rectifier.s_1.v_gs = timeseries(ans.s_r_1.Data(:, 1), ans.s_r_1.Time);
rectifier.s_1.v_ds = timeseries(ans.s_r_1.Data(:, 2), ans.s_r_1.Time);
rectifier.s_1.i = timeseries(ans.s_r_1.Data(:, 3), ans.s_r_1.Time);
rectifier.s_2.v_gs = timeseries(ans.s_r_2.Data(:, 1), ans.s_r_1.Time);
rectifier.s_2.v_ds = timeseries(ans.s_r_2.Data(:, 2), ans.s_r_1.Time);
rectifier.s_2.i = timeseries(ans.s_r_2.Data(:, 3), ans.s_r_1.Time);
rectifier.s_3.v_gs = timeseries(ans.s_r_3.Data(:, 1), ans.s_r_1.Time);
rectifier.s_3.v_ds = timeseries(ans.s_r_3.Data(:, 2), ans.s_r_1.Time);
rectifier.s_3.i = timeseries(ans.s_r_3.Data(:, 3), ans.s_r_1.Time);

rectifier.output.v = timeseries(ans.rectifier_output.Data(:, 1), ans.rectifier_output.Time);
rectifier.output.i = timeseries(ans.rectifier_output.Data(:, 2), ans.rectifier_output.Time);

% 
clear ans;

%% 

save('results/alternator.mat', 'alternator', '-v7.3');
save('results/rectifier.mat', 'rectifier', '-v7.3');

%% Finaliza modelo no Simulink

close_system('AutomotiveEMS.slx');

%% Remoção do diretório principal e seus subdiretórios dos caminhos de 
%  busca do MATLAB

rmpath(root);
