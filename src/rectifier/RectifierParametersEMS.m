%% Parâmetros relativos ao retificador utilizado no Sistema de Gerenciamento
%  de Energia (EMS)

%% Registro de variáveis existente no workspace antes da execução do script

% Registra variáveis existentes no workspace neste instante
if (exist('previousVars', 'var'))
    previousVars{:, end + 1} = who;
else
    previousVars{1} = who;
end

%% Executa script de determinação de parâmetros de acordo com retificador
%  selecionado

switch rectifierCase
    case 'Rivas2004'
        RectifierParametersRivas;
    case 'Sarafianos2015'
        RectifierParametersSarafianos;
    case 'Rodrigues2019'
        disp('Indisponivel');
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
