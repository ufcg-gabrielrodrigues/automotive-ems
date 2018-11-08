%% Parâmetros relativos ao retificador utilizado no Sistema de Gerenciamento
%  de Energia (EMS)

%% Registro de variáveis existente no workspace antes da execução do script

% Registra variáveis existentes no workspace neste instante
if (exist('previousVars', 'var'))
    previousVars{:, end + 1} = who;
else
    previousVars{1} = who;
end

%% Parâmetros de circuito

rectifier.s.v_gs_on = 10;           % Tensão gate-source para acionamento dos MOSFETs do circuito retificador [V]

%% Parâmetros de controle

rectifier.control.pwm.f_s = 10e+3;  % Frequência de chaveamento dos PWMs de controle do circuito retificador [Hz]

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
