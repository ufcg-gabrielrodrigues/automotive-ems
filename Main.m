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

%% Motor a combustão interna

% 

%% Alternador

% Determina realização ou não de nova iteração para determinação de
% parâmetros
alternatorFittingFlag = true;

% Escolha do alternador a ser utilizado
alternatorCase = 'Sarafianos2015';

% Executa script que determina parâmetros do alternador
AlternatorParametersEMS;

% 
k_v_str = regexprep(func2str(alternator.k_v), '@\(.+?\)', '');
replaceFileExpression('models/+SimscapeCustomBlocks/+Controllers/load_matching_smr_controller.ssc', ...
    'k_v == 0', ['k_v == ' k_v_str]);

% 
e_a_str = regexprep(func2str(alternator.stator.input.e.function), '@\(.+?\)', '');
replaceFileExpression('models/+SimscapeCustomBlocks/+Alternator/back_emf.ssc', ...
    'e == 0', ['e == ' e_a_str]);

% 
l_a_str = regexprep(func2str(alternator.stator.l.function), '@\(.+?\)', '');
replaceFileExpression('models/+SimscapeCustomBlocks/+Alternator/stator_inductance.ssc', ...
    'l == 1e-6', ['l == ' l_a_str]);

%% Retificador

% 
RectifierParametersEMS;

%% Curvas P x u

MPPTCurvesRoutine;

%% 

% 
replaceFileExpression('models/+SimscapeCustomBlocks/+Controllers/load_matching_smr_controller.ssc', ...
    ['k_v == ' k_v_str], 'k_v == 0');

% 
replaceFileExpression('models/+SimscapeCustomBlocks/+Alternator/back_emf.ssc', ...
    ['e == ' e_a_str], 'e == 0');

%
replaceFileExpression('models/+SimscapeCustomBlocks/+Alternator/stator_inductance.ssc', ...
    ['l == ' l_a_str], 'l == 1e-6');


%% Remoção do diretório principal e seus subdiretórios dos caminhos de 
%  busca do MATLAB

rmpath(root);
