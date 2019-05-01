%% Análise gráfica para determinação de parâmetros implícitos relativos
%  ao alternador utilizado no artigo de artigo de Sarafianos (2015):
%  "Characterisation and Modelling of Automotive Lundell Alternators"

%% Registro de variáveis existente no workspace antes da execução do script

% Registra variáveis existentes no workspace neste instante
if (exist('previousVars', 'var'))
    previousVars{end + 1} = who;
else
    previousVars{1} = who;
end

%% Identificação dos caminhos de arquivos com dados extraídos dos gráficos

% Caminho absoluto da pasta atual
[currentFolder, ~, ~] = fileparts(mfilename('fullpath'));

% Caminho absoluto da pasta contendo gráficos e dados
figuresFolder = strcat(currentFolder, '/figures');

% Caminhos de subpasta e arquivos referentes aos dados de perdas por atrito
% e enrolamento de acordo com a velocidade do rotor (medições e curvas
% ajustadas)
frictionWindageLossesFolder = strcat(figuresFolder, ...
    '/F&W losses [W] x Speed [rpm]');

frictionWindageLossesMeas = strcat(frictionWindageLossesFolder, ...
    '/F&W losses [W] x Speed [rpm] (Measured points).csv');

frictionWindageLossesFit = strcat(frictionWindageLossesFolder, ...
    '/F&W losses [W] x Speed [rpm] (Curve fitted waveform).csv');

% Caminhos de subpasta e arquivos referentes aos dados de tensão de
% circuito aberto de armadura quando o alternador é submetido a diferentes
% correntes de excitação e velocidades de rotor (medições e curvas
% ajustadas)
openCircuitVoltageFolder = strcat(figuresFolder, ...
    '/OC line-to-line voltage [V] x Field current [A]');

openCircuitVoltage2000rpmMeas = strcat(openCircuitVoltageFolder, ...
    '/OC line-to-line voltage [V] x Field current [A] (2000 rpm - Measured points).csv');
openCircuitVoltage4000rpmMeas = strcat(openCircuitVoltageFolder, ...
    '/OC line-to-line voltage [V] x Field current [A] (4000 rpm - Measured points).csv');
openCircuitVoltage6000rpmMeas = strcat(openCircuitVoltageFolder, ...
    '/OC line-to-line voltage [V] x Field current [A] (6000 rpm - Measured points).csv');
openCircuitVoltage8000rpmMeas = strcat(openCircuitVoltageFolder, ...
    '/OC line-to-line voltage [V] x Field current [A] (8000 rpm - Measured points).csv');

openCircuitVoltage2000rpmFit = strcat(openCircuitVoltageFolder, ...
    '/OC line-to-line voltage [V] x Field current [A] (2000 rpm - Curve fitted waveform).csv');
openCircuitVoltage4000rpmFit = strcat(openCircuitVoltageFolder, ...
    '/OC line-to-line voltage [V] x Field current [A] (4000 rpm - Curve fitted waveform).csv');
openCircuitVoltage6000rpmFit = strcat(openCircuitVoltageFolder, ...
    '/OC line-to-line voltage [V] x Field current [A] (6000 rpm - Curve fitted waveform).csv');
openCircuitVoltage8000rpmFit = strcat(openCircuitVoltageFolder, ...
    '/OC line-to-line voltage [V] x Field current [A] (8000 rpm - Curve fitted waveform).csv');

% Caminhos de subpasta e arquivos referentes aos dados de perdas no ferro
% quando o alternador é submetido a diferentes correntes de excitação e
% velocidades de rotor (medições e curvas ajustadas)
ironLossFolder = strcat(figuresFolder, '/Iron loss [W] x Field current [A]');

ironLoss2000rpmMeas = strcat(ironLossFolder, ...
    '/Iron loss [W] x Field current [A] (2000 rpm - Measured points).csv');
ironLoss4000rpmMeas = strcat(ironLossFolder, ...
    '/Iron loss [W] x Field current [A] (4000 rpm - Measured points).csv');
ironLoss6000rpmMeas = strcat(ironLossFolder, ...
    '/Iron loss [W] x Field current [A] (6000 rpm - Measured points).csv');
ironLoss8000rpmMeas = strcat(ironLossFolder, ...
    '/Iron loss [W] x Field current [A] (8000 rpm - Measured points).csv');

ironLoss2000rpmFit = strcat(ironLossFolder, ...
    '/Iron loss [W] x Field current [A] (2000 rpm - Curve fitted waveform).csv');
ironLoss4000rpmFit = strcat(ironLossFolder, ...
    '/Iron loss [W] x Field current [A] (4000 rpm - Curve fitted waveform).csv');
ironLoss6000rpmFit = strcat(ironLossFolder, ...
    '/Iron loss [W] x Field current [A] (6000 rpm - Curve fitted waveform).csv');
ironLoss8000rpmFit = strcat(ironLossFolder, ...
    '/Iron loss [W] x Field current [A] (8000 rpm - Curve fitted waveform).csv');

% Caminhos de subpasta e arquivos referentes aos dados de indutância de
% armadura quando o alternador é submetido a diferentes correntes de
% excitação (medições e curvas ajustadas)
inductanceFolder = strcat(figuresFolder, '/Inductance [H] x Field Current [A]');

inductanceMeas = strcat(inductanceFolder, ...
    '/Inductance [H] x Field Current [A] (Measured points).csv');

inductanceFit = strcat(inductanceFolder, ...
    '/Inductance [H] x Field Current [A] (Curve fitted waveform).csv');

%% Leitura dos arquivos contendo dados e registro em séries temporais

% Leitura dos dados de perdas por atrito e enrolamento de acordo com a
% velocidade do rotor (medições e curvas ajustadas)
frictionWindageLossesMeas = csvToTimeSeries(frictionWindageLossesMeas);

frictionWindageLossesFit = csvToTimeSeries(frictionWindageLossesFit);

% Leitura dos dados de tensão de circuito aberto de armadura quando o
% alternador é submetido a diferentes correntes de excitação e velocidades
% de rotor (medições e curvas ajustadas)
openCircuitVoltage2000rpmMeas = csvToTimeSeries(openCircuitVoltage2000rpmMeas);
openCircuitVoltage4000rpmMeas = csvToTimeSeries(openCircuitVoltage4000rpmMeas);
openCircuitVoltage6000rpmMeas = csvToTimeSeries(openCircuitVoltage6000rpmMeas);
openCircuitVoltage8000rpmMeas = csvToTimeSeries(openCircuitVoltage8000rpmMeas);

openCircuitVoltage2000rpmFit = csvToTimeSeries(openCircuitVoltage2000rpmFit);
openCircuitVoltage4000rpmFit = csvToTimeSeries(openCircuitVoltage4000rpmFit);
openCircuitVoltage6000rpmFit = csvToTimeSeries(openCircuitVoltage6000rpmFit);
openCircuitVoltage8000rpmFit = csvToTimeSeries(openCircuitVoltage8000rpmFit);

% Leitura dos dados de perdas no ferro quando o alternador é submetido a
% diferentes correntes de excitação e velocidades de rotor (medições e
% curvas ajustadas)
ironLoss2000rpmMeas = csvToTimeSeries(ironLoss2000rpmMeas);
ironLoss4000rpmMeas = csvToTimeSeries(ironLoss4000rpmMeas);
ironLoss6000rpmMeas = csvToTimeSeries(ironLoss6000rpmMeas);
ironLoss8000rpmMeas = csvToTimeSeries(ironLoss8000rpmMeas);

ironLoss2000rpmFit = csvToTimeSeries(ironLoss2000rpmFit);
ironLoss4000rpmFit = csvToTimeSeries(ironLoss4000rpmFit);
ironLoss6000rpmFit = csvToTimeSeries(ironLoss6000rpmFit);
ironLoss8000rpmFit = csvToTimeSeries(ironLoss8000rpmFit);

% Leitura dos dados de indutância de armadura quando o alternador é
% submetido a diferentes correntes de excitação (medições e curvas
% ajustadas)
inductanceMeas = csvToTimeSeries(inductanceMeas);

inductanceFit = csvToTimeSeries(inductanceFit);

%% 

if (alternatorParamPlotFlag)
    figure(1)
    plot(frictionWindageLossesMeas.time, frictionWindageLossesMeas.data, 'bo'); hold on;
    xlabel('$n_{r}\,[rpm]$');
    ylabel('$p_{f\&w}\,[W]$');
    grid on;
    
    figure(2)
    plot(openCircuitVoltage2000rpmMeas.time, openCircuitVoltage2000rpmMeas.data, 'ro'); hold on;
    plot(openCircuitVoltage4000rpmMeas.time, openCircuitVoltage4000rpmMeas.data, 'go'); hold on;
    plot(openCircuitVoltage6000rpmMeas.time, openCircuitVoltage6000rpmMeas.data, 'bo'); hold on;
    plot(openCircuitVoltage8000rpmMeas.time, openCircuitVoltage8000rpmMeas.data, 'ko'); hold on;
    xlabel('$i_{f}\,[A]$');
    ylabel('$e_{ll}\,[V]$');
    grid on;
    
    figure(3)
    plot(ironLoss2000rpmMeas.time, ironLoss2000rpmMeas.data, 'ro'); hold on;
    plot(ironLoss4000rpmMeas.time, ironLoss4000rpmMeas.data, 'go'); hold on;
    plot(ironLoss6000rpmMeas.time, ironLoss6000rpmMeas.data, 'bo'); hold on;
    plot(ironLoss8000rpmMeas.time, ironLoss8000rpmMeas.data, 'ko'); hold on;
    xlabel('$i_{f}\,[A]$');
    ylabel('$p_{i}\,[W]$');
    grid on;
    
    figure(4)
    plot(inductanceMeas.time, inductanceMeas.data, 'bo'); hold on;
    xlabel('$i_{f}\,[A]$');
    ylabel('$l_{s}\,[H]$');
    grid on;
end

%% Registro dos dados extraídos dos gráficos em arquivos .MAT

% Registro dos dados de perdas por atrito e enrolamento de acordo com a
% velocidade do rotor (medições e curvas ajustadas)
registerPath = strcat(frictionWindageLossesFolder, '/frictionWindageLossesMeas.mat');
save(registerPath, 'frictionWindageLossesMeas');

registerPath = strcat(frictionWindageLossesFolder, '/frictionWindageLossesFit.mat');
save(registerPath, 'frictionWindageLossesFit');

% Registro dos dados de tensão de circuito aberto de armadura quando o
% alternador é submetido a diferentes correntes de excitação e velocidades
% de rotor (medições e curvas ajustadas)
registerPath = strcat(openCircuitVoltageFolder, '/openCircuitVoltageMeas.mat');
save(registerPath, 'openCircuitVoltage2000rpmMeas', ...
    'openCircuitVoltage4000rpmMeas', 'openCircuitVoltage6000rpmMeas', ...
    'openCircuitVoltage8000rpmMeas');

registerPath = strcat(openCircuitVoltageFolder, '/openCircuitVoltageFit.mat');
save(registerPath, 'openCircuitVoltage2000rpmFit', ...
    'openCircuitVoltage4000rpmFit', 'openCircuitVoltage6000rpmFit', ...
    'openCircuitVoltage8000rpmFit');

% Registro dos dados de perdas no ferro quando o alternador é submetido a
% diferentes correntes de excitação e velocidades de rotor (medições e
% curvas ajustadas)
registerPath = strcat(ironLossFolder, '/ironLossMeas.mat');
save(registerPath, 'ironLoss2000rpmMeas', 'ironLoss4000rpmMeas', ...
    'ironLoss6000rpmMeas', 'ironLoss8000rpmMeas');

registerPath = strcat(ironLossFolder, '/ironLossFit.mat');
save(registerPath, 'ironLoss2000rpmFit', 'ironLoss4000rpmFit', ...
    'ironLoss6000rpmFit', 'ironLoss8000rpmFit');

% Registro dos dados de indutância de armadura quando o alternador é
% submetido a diferentes correntes de excitação (medições e curvas
% ajustadas)
registerPath = strcat(inductanceFolder, '/inductanceMeas.mat');
save(registerPath, 'inductanceMeas');

registerPath = strcat(inductanceFolder, '/inductanceFit.mat');
save(registerPath, 'inductanceFit');

%% Preparação de dados medidos para etapa de ajustes de curvas

frictionWindageLossesMeas = scaleMixTimeSeries(frictionWindageLossesMeas, 1.0);

openCircuitVoltageMeas = scaleMixTimeSeries(openCircuitVoltage2000rpmMeas, 1/2000, ...
    openCircuitVoltage4000rpmMeas, 1/4000, openCircuitVoltage6000rpmMeas, 1/6000, ...
    openCircuitVoltage8000rpmMeas, 1/8000);

ironLossMeas = scaleMixTimeSeries(ironLoss2000rpmMeas, 1/2000, ...
    ironLoss4000rpmMeas, 1/4000, ironLoss6000rpmMeas, 1/6000, ...
    ironLoss8000rpmMeas, 1/8000);

inductanceMeas = scaleMixTimeSeries(inductanceMeas, 1.0);

%% Ajuste de curvas segundo método apropriado

% Tipos de ajuste de curva
quadraticFit = 'poly2';     % Ajuste por equação de segunda ordem
cubicFit = 'poly3';         % Ajuste por equação de terceira ordem
fourthOrderFit = 'poly4';   % Ajuste por equação de quarta ordem

% Realização do ajuste de curvas segundo tipo escolhido para cada caso
frictionWindageLossesFit = fitTimeSeries(frictionWindageLossesMeas, quadraticFit);

openCircuitVoltageFit = fitTimeSeries(openCircuitVoltageMeas, fourthOrderFit);

ironLossFit = fitTimeSeries(ironLossMeas, fourthOrderFit);

inductanceFit = fitTimeSeries(inductanceMeas, cubicFit);

%% Utilização dos parâmetros encontrados via ajuste de curvas para criar
%  manipuladores de função

frictionWindageLosses = vectorToPolyFunctionHandle(coeffvalues(frictionWindageLossesFit), 'omega_r');

openCircuitVoltage = vectorToPolyFunctionHandle(coeffvalues(openCircuitVoltageFit), 'i_f');

ironLoss = vectorToPolyFunctionHandle(coeffvalues(ironLossFit), 'i_f');

inductance = vectorToPolyFunctionHandle(coeffvalues(inductanceFit), 'i_f');

%% Registro das características implícitas em arquivo .MAT

% Concatenação do caminho absoluto da pasta atual com nome do arquivo
implicitAlternatorCharacteristics = strcat(currentFolder, ...
    '/implicitAlternatorCharacteristics.mat');

% Registro das características no destino previamente especificado
save(implicitAlternatorCharacteristics, 'frictionWindageLosses', ...
    'openCircuitVoltage', 'ironLoss', 'inductance');

%% 

if (alternatorParamPlotFlag)
    figure(1)
    fplot(frictionWindageLosses, [2000 8000], 'b-'); hold off;
    
    figure(2)
    fplot(@(i_f) 2000*openCircuitVoltage(i_f), [0 5], 'r-'); hold on;
    fplot(@(i_f) 4000*openCircuitVoltage(i_f), [0 5], 'g-'); hold on;
    fplot(@(i_f) 6000*openCircuitVoltage(i_f), [0 5], 'b-'); hold on;
    fplot(@(i_f) 8000*openCircuitVoltage(i_f), [0 5], 'k-'); hold off;
    
    figure(3)
    fplot(@(i_f) 2000*ironLoss(i_f), [0 5], 'r-'); hold on;
    fplot(@(i_f) 4000*ironLoss(i_f), [0 5], 'g-'); hold on;
    fplot(@(i_f) 6000*ironLoss(i_f), [0 5], 'b-'); hold on;
    fplot(@(i_f) 8000*ironLoss(i_f), [0 5], 'k-'); hold off;
    
    figure(4)
    fplot(inductance, [0 5], 'b-'); hold off;
    
    saveFigure(figure(1), 'results/LundellAlternator/friction-windage-loss', 'fig');
    saveFigure(figure(1), 'results/LundellAlternator/friction-windage-loss', 'png');
    
    saveFigure(figure(2), 'results/LundellAlternator/open-circuit-voltage', 'fig');
    saveFigure(figure(2), 'results/LundellAlternator/open-circuit-voltage', 'png');
    
    saveFigure(figure(3), 'results/LundellAlternator/iron-loss', 'fig');
    saveFigure(figure(3), 'results/LundellAlternator/iron-loss', 'png');
    
    saveFigure(figure(4), 'results/LundellAlternator/inductance', 'fig');
    saveFigure(figure(4), 'results/LundellAlternator/inductance', 'png');
end

%% Exclusão das variáveis excedentes

% Inicializa variáveis para que apareçam no workspace
currentVars = [];
newVars = [];

% Identifica variáveis existentes no workspace neste instante
currentVars = who;

% Pela diferença, determina variáveis criadas pelo script e, se necessário,
% apaga registro atual de variáveis
newVars = setdiff(currentVars, previousVars{end});

% Exclui último registro de variáveis do workspace
if (length(previousVars) > 1)
    previousVars(end) = [];
end

% Elimina exceções da lista de variáveis a serem deletadas
newVars(multStrCmp(newVars, {'implicitAlternatorCharacteristics'})) = [];

% Realiza limpeza
clear(newVars{:});
