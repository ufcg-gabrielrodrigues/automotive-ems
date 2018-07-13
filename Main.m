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

open_system('AutomotiveEMS.slx');

%% Motor a combustão interna

ice.speed = 2000;       % Velocidade de rotação [rpm]

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
    sym(alternator.e_a), 'FunctionName', 'backEMF', 'Outputs', {'e_a'});

% 
matlabFunctionBlock('AutomotiveEMS/Alternator/Stator Inductance (abc reference)/statorInductance', ...
    sym(alternator.stator.l), 'FunctionName', 'statorInductance', 'Outputs', {'l_s'});

%% Retificador

% 
RectifierParametersEMS;

%% 

iceToAltRotRatio = 2.5;

%% 

T_s = 1e-6;

%% Remoção do diretório principal e seus subdiretórios dos caminhos de 
%  busca do MATLAB

rmpath(root);
