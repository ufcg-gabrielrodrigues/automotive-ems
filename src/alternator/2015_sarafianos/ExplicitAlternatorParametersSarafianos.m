%% Características explícitas relativos ao alternador utilizado no artigo de
%  Sarafianos (2015): "Characterisation and Modelling of Automotive Lundell 
%  Alternators"

%% Registro de variáveis existente no workspace antes da execução do script

% Registra variáveis existentes no workspace neste instante
if (exist('previousVars', 'var'))
    previousVars{:, end + 1} = who;
else
    previousVars{1} = who;
end

%% Características explícitas do alternador utilizado

n = 3;          % Número de fases
p = 8;          % Número de pares de polos por fase
slots = 48;     % Número de ranhuras no estator
r_a = 0.03;     % Resistência de estator (por fase) a 20 oC
r_f = 1.90;     % Resistência de rotor a 20 oC
v_o = 13.5;     % Tensão de saída (CC)
i_max_2k = 120; % Corrente máxima de saída a 2000 rpm e 25 oC
i_max_4k = 170; % Corrente máxima de saída a 4000 rpm e 25 oC
i_max_6k = 180; % Corrente máxima de saída a 6000 rpm e 25 oC

%% Registro das características explícitas em arquivo .MAT

% Caminho absoluto da pasta atual onde serão registradas as características
[currentFolder, ~, ~] = fileparts(mfilename('fullpath'));

% Concatenação do caminho absoluto da pasta atual com nome do arquivo
explicitAlternatorCharacteristics = strcat(currentFolder, ...
    '/explicitAlternatorCharacteristics.mat');

% Registro das características no destino previamente especificado
save(explicitAlternatorCharacteristics, 'n', 'p', 'slots', 'r_a', 'r_f', 'v_o', ...
    'i_max_2k', 'i_max_4k', 'i_max_6k');

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
newVars(multStrCmp(newVars, {'explicitAlternatorCharacteristics'})) = [];

% Realiza limpeza
clear(newVars{:});
