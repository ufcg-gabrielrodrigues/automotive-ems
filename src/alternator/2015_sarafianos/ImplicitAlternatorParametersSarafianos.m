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
    xlabel('$n_{r}\,[\textrm{rpm}]$');
    ylabel('$p_{f\&w}\,[\textrm{W}]$');
    grid on;
    
    figure(2)
    plot(openCircuitVoltage2000rpmMeas.time, openCircuitVoltage2000rpmMeas.data, 'ro'); hold on;
    plot(openCircuitVoltage4000rpmMeas.time, openCircuitVoltage4000rpmMeas.data, 'go'); hold on;
    plot(openCircuitVoltage6000rpmMeas.time, openCircuitVoltage6000rpmMeas.data, 'bo'); hold on;
    plot(openCircuitVoltage8000rpmMeas.time, openCircuitVoltage8000rpmMeas.data, 'ko'); hold on;
    xlabel('$i_{f}\,[\textrm{A}]$');
    ylabel('$e_{s,\textrm{rms}}^{ll}\,[\textrm{V}]$');
    legend('$n_{r} = 2000\,\textrm{rpm}$', '$n_{r} = 4000\,\textrm{rpm}$', ...
           '$n_{r} = 6000\,\textrm{rpm}$', '$n_{r} = 8000\,\textrm{rpm}$', ...
           'Location', 'NorthWest');
    grid on;
    
    figure(3)
    plot(ironLoss2000rpmMeas.time, ironLoss2000rpmMeas.data, 'ro'); hold on;
    plot(ironLoss4000rpmMeas.time, ironLoss4000rpmMeas.data, 'go'); hold on;
    plot(ironLoss6000rpmMeas.time, ironLoss6000rpmMeas.data, 'bo'); hold on;
    plot(ironLoss8000rpmMeas.time, ironLoss8000rpmMeas.data, 'ko'); hold on;
    xlabel('$i_{f}\,[\textrm{A}]$');
    ylabel('$p_{i}\,[\textrm{W}]$');
    legend('$n_{r} = 2000\,\textrm{rpm}$', '$n_{r} = 4000\,\textrm{rpm}$', ...
           '$n_{r} = 6000\,\textrm{rpm}$', '$n_{r} = 8000\,\textrm{rpm}$', ...
           'Location', 'NorthWest');
    grid on;
    
    figure(4)
    plot(inductanceMeas.time, inductanceMeas.data, 'bo'); hold on;
    xlabel('$i_{f}\,[\textrm{A}]$');
    ylabel('$l_{s}\,[\textrm{H}]$');
    grid on;
    
    saveFigure(figure(1), 'results/LundellAlternator/friction-windage-loss-exp', 'fig');
    saveFigure(figure(1), 'results/LundellAlternator/friction-windage-loss-exp', 'eps');
    
    saveFigure(figure(2), 'results/LundellAlternator/open-circuit-voltage-exp', 'fig');
    saveFigure(figure(2), 'results/LundellAlternator/open-circuit-voltage-exp', 'eps');
    
    saveFigure(figure(3), 'results/LundellAlternator/iron-loss-exp', 'fig');
    saveFigure(figure(3), 'results/LundellAlternator/iron-loss-exp', 'eps');
    
    saveFigure(figure(4), 'results/LundellAlternator/inductance-exp', 'fig');
    saveFigure(figure(4), 'results/LundellAlternator/inductance-exp', 'eps');
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

%% Ajuste de curvas segundo método apropriado

% 
frictionWindageLossesFit = fitTimeSeries(frictionWindageLossesMeas, 'poly2');
frictionWindageLosses = vectorToPolyFunctionHandle(frictionWindageLossesFit, 'omega_r');

% 
openCircuitVoltageFit(1, :) = fitTimeSeries(openCircuitVoltage2000rpmMeas, 'sigmoid');
openCircuitVoltageFit(2, :) = fitTimeSeries(openCircuitVoltage4000rpmMeas, 'sigmoid');
openCircuitVoltageFit(3, :) = fitTimeSeries(openCircuitVoltage6000rpmMeas, 'sigmoid');
openCircuitVoltageFit(4, :) = fitTimeSeries(openCircuitVoltage8000rpmMeas, 'sigmoid');

openCircuitVoltage{1} = vectorToSigmFunctionHandle(openCircuitVoltageFit(1, :), 'i_f');
openCircuitVoltage{2} = vectorToSigmFunctionHandle(openCircuitVoltageFit(2, :), 'i_f');
openCircuitVoltage{3} = vectorToSigmFunctionHandle(openCircuitVoltageFit(3, :), 'i_f');
openCircuitVoltage{4} = vectorToSigmFunctionHandle(openCircuitVoltageFit(4, :), 'i_f');

i_f = 0.00:0.25:5.00;
n_r = 2000:2000:8000;
e_ll = zeros(length(n_r), length(i_f));

e_ll(1, :) = openCircuitVoltage{1}(i_f);
e_ll(2, :) = openCircuitVoltage{2}(i_f);
e_ll(3, :) = openCircuitVoltage{3}(i_f);
e_ll(4, :) = openCircuitVoltage{4}(i_f);

openCircuitVoltageFun = 'n_r.*i_f.*(a+(b-a)./(1+10.^((c-i_f)*d)))';

openCircuitVoltage = fitToFunction(customSurfaceFit(i_f, n_r, e_ll, openCircuitVoltageFun, 'i_f', 'n_r', 'e_ll', [0.292136752238811 0.165676040245502 0.164950588884314 0.906364323265215]));

% 
ironLossFit(1, :) = fitTimeSeries(ironLoss2000rpmMeas, 'poly7');
ironLossFit(2, :) = fitTimeSeries(ironLoss4000rpmMeas, 'poly7');
ironLossFit(3, :) = fitTimeSeries(ironLoss6000rpmMeas, 'poly7');
ironLossFit(4, :) = fitTimeSeries(ironLoss8000rpmMeas, 'poly7');

ironLoss{1} = vectorToPolyFunctionHandle(ironLossFit(1, :), 'i_f');
ironLoss{2} = vectorToPolyFunctionHandle(ironLossFit(2, :), 'i_f');
ironLoss{3} = vectorToPolyFunctionHandle(ironLossFit(3, :), 'i_f');
ironLoss{4} = vectorToPolyFunctionHandle(ironLossFit(4, :), 'i_f');

i_f = 0.00:0.25:5.00;
n_r = 2000:2000:8000;
p_i = zeros(length(n_r), length(i_f));

p_i(1, :) = ironLoss{1}(i_f);
p_i(2, :) = ironLoss{2}(i_f);
p_i(3, :) = ironLoss{3}(i_f);
p_i(4, :) = ironLoss{4}(i_f);

ironLossFun = '(k_n0 + k_n1.*n_r).*(k_i0 + k_i1.*(i_f) + k_i2.*(i_f.^2) + k_i3.*(i_f.^3) + k_i4.*(i_f.^4) + k_i5.*(i_f.^5) + k_i6.*(i_f.^6) + k_i7.*(i_f.^7))';

ironLoss = fitToFunction(customSurfaceFit(i_f, n_r, p_i, ironLossFun, 'i_f', 'n_r', 'p_i', [0.720855670816931 0.361022049194661 0.620278427071085 0.811150885100285 0.0192574774141414 0.0838735082828999 0.97480166718489 0.651349532415353 0.0186127747263861 0.231237816164352]));

% 
inductanceFit = fitTimeSeries(inductanceMeas, 'poly3');
inductance = vectorToPolyFunctionHandle(inductanceFit, 'i_f');

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
    x = 2000:1:8000;
    
    y = frictionWindageLosses(x);
    plot(x, y, 'b-'); hold off;
    
    figure(2)
    x = 0:1e-3:5;
    
    y = openCircuitVoltage(x, 2000);
    plot(x, y, 'r-', 'HandleVisibility', 'off'); hold on;
    
    y = openCircuitVoltage(x, 4000);
    plot(x, y, 'g-', 'HandleVisibility', 'off'); hold on;
    
    y = openCircuitVoltage(x, 6000);
    plot(x, y, 'b-', 'HandleVisibility', 'off'); hold on;
    
    y = openCircuitVoltage(x, 8000);
    plot(x, y, 'k-', 'HandleVisibility', 'off'); hold off;
    
    figure(3)
    x = 0:1e-3:5;
    
    y = ironLoss(x, 2000);
    plot(x, y, 'r-', 'HandleVisibility', 'off'); hold on;
    
    y = ironLoss(x, 4000);
    plot(x, y, 'g-', 'HandleVisibility', 'off'); hold on;
    
    y = ironLoss(x, 6000);
    plot(x, y, 'b-', 'HandleVisibility', 'off'); hold on;
    
    y = ironLoss(x, 8000);
    plot(x, y, 'k-', 'HandleVisibility', 'off'); hold off;
    
    figure(4)
    legend('off');
    x = 0:1e-3:5;
    
    y = inductance(x);
    plot(x, y, 'b-', 'HandleVisibility', 'off'); hold off;
    
    saveFigure(figure(1), 'results/LundellAlternator/friction-windage-loss-fit', 'fig');
    saveFigure(figure(1), 'results/LundellAlternator/friction-windage-loss-fit', 'eps');
    
    saveFigure(figure(2), 'results/LundellAlternator/open-circuit-voltage-fit', 'fig');
    saveFigure(figure(2), 'results/LundellAlternator/open-circuit-voltage-fit', 'eps');
    
    saveFigure(figure(3), 'results/LundellAlternator/iron-loss-fit', 'fig');
    saveFigure(figure(3), 'results/LundellAlternator/iron-loss-fit', 'eps');
    
    saveFigure(figure(4), 'results/LundellAlternator/inductance-fit', 'fig');
    saveFigure(figure(4), 'results/LundellAlternator/inductance-fit', 'eps');
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
