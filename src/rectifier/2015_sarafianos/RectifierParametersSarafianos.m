%% Características relativos ao retificador utilizado no artigo de
%  Sarafianos (2015): "Characterisation and Modelling of Automotive Lundell
%  Alternators"

%% Registro de variáveis existente no workspace antes da execução do script

% Registra variáveis existentes no workspace neste instante
if (exist('previousVars', 'var'))
    previousVars{:, end + 1} = who;
else
    previousVars{1} = who;
end

%% Realiza nova iteração de análise das características do retificador ou
%  apenas carrega os resultados de uma iteração passada

if (rectifierCalcParamFlag)
    %% Registro de parâmetros definidos em estrutura que representa o retificador
    
    % Diodos
    rectifier.d.v_on = 0.8;             % Tensão de polarização [V]
    rectifier.d.r_on = 10e-3;           % Resistência de condução [Ohm]
    
    % MOSFETs
    rectifier.s.v_gs_on = 10;           % Tensão gate-source para acionamento dos MOSFETs do circuito retificador [V]
    
    % Filtro passivo
    rectifier.filter.c = 47e-3;         % Capacitância de filtro [F]
    
    % Parâmetros de controle
    rectifier.control.pwm.f_s = 1.0e+4; % Frequência de chaveamento dos PWMs de controle do circuito retificador [Hz]
    
    %% Registro da estrutura que representa o retificador em arquivos .MAT
    
    % Caminho absoluto da pasta atual onde será registrada a estrutura
    [currentFolder, ~, ~] = fileparts(mfilename('fullpath'));
    
    % Concatenação do caminho absoluto das pastas atual e acima com nome do
    % arquivo
    registerCurrentPath = strcat(currentFolder, '/rectifier.mat');
    registerAbovePath = strcat(currentFolder, '/../rectifier.mat');
    
    % Registro da estrutura nos destinos previamente especificados
    save(registerCurrentPath, 'rectifier');
    save(registerAbovePath, 'rectifier');
else
    %% Carrega estrutura que representa o retificador obtida em iterações
    %  anteriores e armazena cópia no diretório acima
    
    % Caminho absoluto da pasta atual onde será registrada a estrutura
    [currentFolder, ~, ~] = fileparts(mfilename('fullpath'));
    
    % Concatenação do caminho absoluto das pastas atual e acima com nome do
    % arquivo
    registerCurrentPath = strcat(currentFolder, '/rectifier.mat');
    registerAbovePath = strcat(currentFolder, '/../rectifier.mat');
    
    % Registro da estrutura nos destinos previamente especificados
    load(registerCurrentPath);
    save(registerAbovePath, 'rectifier');
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
newVars(multStrCmp(newVars, {'rectifier'})) = [];

% Realiza limpeza
clear(newVars{:});
