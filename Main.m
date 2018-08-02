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

% 
e_a_str = regexprep(func2str(alternator.stator.input.e.function), '@\(.+?\)', '');
replaceFileExpression('+SimscapeCustomBlocks/+Alternator/back_emf.ssc', ...
    'e == 0', ['e == ' e_a_str]);

% 
l_a_str = regexprep(func2str(alternator.stator.l.function), '@\(.+?\)', '');
replaceFileExpression('+SimscapeCustomBlocks/+Alternator/stator_inductance.ssc', ...
    'l == 1e-6', ['l == ' l_a_str]);

%% Retificador

% 
RectifierParametersEMS;

%% 

iceToAltRotRatio = 2.5;

%% 

t = 1e+1;           % [s]

%% 

simulationParameters.StopTime   = num2str(t);

%% 

sim('AutomotiveEMS', simulationParameters);

%% 

% Motor a combustão interna
ice.n = ans.n_ice;

% Alternador
alternator.rotor.control.q = ans.q_s_f;
alternator.rotor.l.i = ans.i_f;

alternator.stator.input.e.value = ans.e_a_abc;
alternator.stator.l.value = ans.l_a_abc;
alternator.stator.output.v = ans.v_a_abc;

% Retificador
rectifier.output.v = timeseries(ans.rectifier_output.Data(:, 1), ans.rectifier_output.Time);
rectifier.output.i = timeseries(ans.rectifier_output.Data(:, 2), ans.rectifier_output.Time);

% 
clear ans;

%% 

save('results/alternator.mat', 'alternator', '-v7.3');
save('results/rectifier.mat', 'rectifier', '-v7.3');

%% Finaliza modelo no Simulink

close_system('AutomotiveEMS.slx');

%% 

% 
replaceFileExpression('+SimscapeCustomBlocks/+Alternator/back_emf.ssc', ...
    ['e == ' e_a_str], 'e == 0');

%
replaceFileExpression('+SimscapeCustomBlocks/+Alternator/stator_inductance.ssc', ...
    ['l == ' l_a_str], 'l == 1e-6');

%% Remoção do diretório principal e seus subdiretórios dos caminhos de 
%  busca do MATLAB

rmpath(root);
