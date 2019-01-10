%% Características relativos ao alternador (Ford 130 A) utilizado no
%  artigo de Rivas (2004): "Performance Improvement of Alternators With
%  Switched-Mode Rectifiers"

%% Registro de variáveis existente no workspace antes da execução do script

% Registra variáveis existentes no workspace neste instante
if (exist('previousVars', 'var'))
    previousVars{:, end + 1} = who;
else
    previousVars{1} = who;
end

%% Realiza nova iteração de análise das características do alternador ou
%  apenas carrega os resultados de uma iteração passada

if (alternatorCalcParamFlag)
    %% Registro de parâmetros definidos em estrutura que representa o alternador
    
    alternator.n = 3;                                           % Número de fases
    alternator.p = 6;                                           % Número de pares de polos por fase
    alternator.stator.connection = y;                           % Tipo de conexão do circuito de estator
    alternator.stator.l.value = 120e-6;                         % Indutância própria por fase do circuito de armadura [H]
    alternator.stator.r.value = 37e-3;                          % Resistência por fase do circuito de armadura [Ohm]
    alternator.k_e.value = 10.716*sqrt(2)/((2*pi*180)*(3.6));	% Constante de acoplamento el�trico [V/((rad/s)*A)]
    
    %% Registro da estrutura que representa o alternador em arquivos .MAT
    
    % Caminho absoluto da pasta atual onde será registrada a estrutura
    [currentFolder, ~, ~] = fileparts(mfilename('fullpath'));
    
    % Concatenação do caminho absoluto das pastas atual e acima com nome do
    % arquivo
    registerCurrentPath = strcat(currentFolder, '/alternator.mat');
    registerAbovePath = strcat(currentFolder, '/../alternator.mat');
    
    % Registro da estrutura nos destinos previamente especificados
    save(registerCurrentPath, 'alternator');
    save(registerAbovePath, 'alternator');
else
    %% Carrega estrutura que representa o alternador obtida em iterações
    %  anteriores e armazena cópia no diretório acima
    
    % Caminho absoluto da pasta atual onde será registrada a estrutura
    [currentFolder, ~, ~] = fileparts(mfilename('fullpath'));
    
    % Concatenação do caminho absoluto das pastas atual e acima com nome do
    % arquivo
    registerCurrentPath = strcat(currentFolder, '/alternator.mat');
    registerAbovePath = strcat(currentFolder, '/../alternator.mat');
    
    % Registro da estrutura nos destinos previamente especificados
    load(registerCurrentPath);
    save(registerAbovePath, 'alternator');
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
newVars(multStrCmp(newVars, {'alternator'})) = [];

% Realiza limpeza
clear(newVars{:});
