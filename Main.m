%% Script principal da simulação do Sistema de Gerenciamento de Energia (EMS)

%% Preparação de ambiente para simulação

clear all;
close all;
clc;

%% Formatação de gráficos

% Formatação
set(groot, 'DefaultTextInterpreter', 'LaTeX');
set(groot, 'DefaultLegendInterpreter', 'LaTeX');
set(groot, 'DefaultTextFontSize', 25);
set(groot, 'DefaultAxesFontSize', 25);
set(groot, 'DefaultLegendFontSize', 25);
set(groot, 'DefaultLineLineWidth', 2.5);

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

% Tipos de conexão do circuito de estator
y = 1;
delta = 2;

% Determina realização ou não de nova iteração para determinação de
% parâmetros
alternatorCalcParamFlag = true;

% Determina realização ou não do traçado dos resultados da determinação de
% parâmetros
alternatorParamPlotFlag = true;

% Determina realização ou não do traçado dos parâmetros definidos por
% função
alternatorFuncParamPlotFlag = true;

% Escolha do alternador a ser utilizado
alternatorCase = 'Sarafianos2015';

% Executa script que determina parâmetros do alternador
AlternatorParametersEMS;

% Atualização de parâmetro: indutância mútua entre armadura e campo
m_f_default = 'm_f == { 0, ''V/((rad/s)*(A))'' };';

if (isfield(alternator.m_f, 'function'))
    m_f_str = regexprep(func2str(alternator.m_f.function), '@\(.+?\)', '');
    m_f_str = strrep(m_f_str, 'i_f', '(i_f*{1,''1/A''})');
    
    if (alternatorFuncParamPlotFlag)
        figure(1)
        x = 0:1e-3:alternator.rotor.i_max;
        y = alternator.m_f.function(x);
        plot(x, y);
        xlabel('$i_{f}\,[A]$');
        ylabel('$m_{f}\,[\frac{V\,s}{A\,rad}]$');
        grid on;
        
        fileName = sprintf('results/LundellAlternator/alternador-indutancia-mutua');
        saveFigure(figure(1), fileName, 'fig');
        saveFigure(figure(1), fileName, 'eps');
        close all;
    end
else
    m_f_str = num2str(alternator.m_f.value);
end

m_f = ['m_f == { ' m_f_str ', ''V/((rad/s)*(A))'' };'];
replaceFileExpression('models/+SimscapeCustomBlocks/+Alternator/back_emf.ssc', ...
    m_f_default, m_f);

% Atualização de parâmetro: indutância própria de estator
l_s_default = 'l == { 1e-6, ''H'' };';

if (isfield(alternator.stator.l, 'function'))
    l_s_str = regexprep(func2str(alternator.stator.l.function), '@\(.+?\)', '');
    l_s_str = strrep(l_s_str, 'i_f', '(i_f*{1,''1/A''})');
    
    if (alternatorFuncParamPlotFlag)
        figure(1)
        x = 0:1e-3:alternator.rotor.i_max;
        y = alternator.stator.l.function(x);
        plot(x, y);
        xlabel('$i_{f}\,[A]$');
        ylabel('$l_{s}\,[H]$');
        grid on;
        
        fileName = sprintf('results/LundellAlternator/alternador-indutancia-estator');
        saveFigure(figure(1), fileName, 'fig');
        saveFigure(figure(1), fileName, 'eps');
        close all;
    end
else
    l_s_str = num2str(alternator.stator.l.value);
end

l_s = ['l == { ' l_s_str ', ''H'' };'];
replaceFileExpression('models/+SimscapeCustomBlocks/+Alternator/stator_inductance.ssc', ...
    l_s_default, l_s);

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

LundellAlternatorRoutine;

%% Redefinição de parâmetros de alternador

% Atualização de parâmetro para valor padrão: indutância mútua entre
% armadura e campo
replaceFileExpression('models/+SimscapeCustomBlocks/+Alternator/back_emf.ssc', ...
    m_f, m_f_default);

% Atualização de parâmetro para valor padrão: indutância própria de estator
replaceFileExpression('models/+SimscapeCustomBlocks/+Alternator/stator_inductance.ssc', ...
    l_s, l_s_default);

%% Remoção do diretório principal e seus subdiretórios dos caminhos de 
%  busca do MATLAB

rmpath(root);
