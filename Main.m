%% Script principal da simulação do Sistema de Gerenciamento de Energia (EMS)

%% Preparação de ambiente para simulação

clear all;
close all;
clc;

%% Formatação de gráficos

% Formatação
set(groot, 'DefaultTextInterpreter', 'LaTeX');
set(groot, 'DefaultLegendInterpreter', 'LaTeX');
set(groot, 'DefaultTextFontSize', 14);
set(groot, 'DefaultAxesFontSize', 14);
set(groot, 'DefaultLegendFontSize', 14);
set(groot, 'DefaultLineLineWidth', 1);

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
alternatorCalcParamFlag = true;

% Escolha do alternador a ser utilizado
alternatorCase = 'Sarafianos2015';

% Executa script que determina parâmetros do alternador
AlternatorParametersEMS;

% Atualização de parâmetro: constante de acoplamento el�trico
k_e_default = 'k_e == { 0, ''V/((rad/s)*(A))'' };';

if (isfield(alternator.k_e, 'function'))
    k_e_str = regexprep(func2str(alternator.k_e.function), '@\(.+?\)', '');
    k_e_str = strrep(k_e_str, 'i_f', '(i_f*{1,''1/A''})');
else
    k_e_str = num2str(alternator.k_e.value);
end

k_e = ['k_e == { ' k_e_str ', ''V/((rad/s)*(A))'' };'];
replaceFileExpression('models/+SimscapeCustomBlocks/+Alternator/back_emf.ssc', ...
    k_e_default, k_e);

% Atualização de parâmetro: indutância própria de estator
l_a_default = 'l == { 1e-6, ''H'' };';

if (isfield(alternator.stator.l, 'function'))
    l_a_str = regexprep(func2str(alternator.stator.l.function), '@\(.+?\)', '');
    l_a_str = strrep(l_a_str, 'i_f', '(i_f*{1,''1/A''})');
else
    l_a_str = num2str(alternator.stator.l.value);
end

l_a = ['l == { ' l_a_str ', ''H'' };'];
replaceFileExpression('models/+SimscapeCustomBlocks/+Alternator/stator_inductance.ssc', ...
    l_a_default, l_a);

% Tipos de conexão do circuito de estator
yStatorConnection = Simulink.Variant('alternator.stator.connection == 1');
deltaStatorConnection = Simulink.Variant('alternator.stator.connection == 2');

%% Retificador

% Determina realização ou não de nova iteração para determinação de
% parâmetros
rectifierCalcParamFlag = true;

% Escolha do retificador a ser utilizado
rectifierCase = 'Sarafianos2015';

% Executa script que determina parâmetros do retificador
RectifierParametersEMS;

%% Rotina de simulação

MPPTCurvesRoutine;

%% Redefinição de parâmetros de alternador

% Atualização de parâmetro para valor padrão: constante de acoplamento
% el�trico
replaceFileExpression('models/+SimscapeCustomBlocks/+Alternator/back_emf.ssc', ...
    k_e, k_e_default);

% Atualização de parâmetro para valor padrão: indutância própria de estator
replaceFileExpression('models/+SimscapeCustomBlocks/+Alternator/stator_inductance.ssc', ...
    l_a, l_a_default);

%% Remoção do diretório principal e seus subdiretórios dos caminhos de 
%  busca do MATLAB

rmpath(root);
