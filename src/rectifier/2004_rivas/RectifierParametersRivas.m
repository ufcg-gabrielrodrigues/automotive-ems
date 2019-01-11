%% Características relativos ao retificador utilizado no artigo de Rivas
%  (2004): "Performance Improvement of Alternators With Switched-Mode
%  Rectifiers"

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
    
    % MOSFETs
    rectifier.s.v_gs_on = 10;           % Tensão gate-source para acionamento dos MOSFETs do circuito retificador [V]
    
    % Parâmetros de controle
    rectifier.control.pwm.f_s = 1.0e+5; % Frequência de chaveamento dos PWMs de controle do circuito retificador [Hz]
    
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
