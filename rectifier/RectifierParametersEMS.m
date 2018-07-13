%% Parâmetros relativos ao retificador utilizado no Sistema de Gerenciamento
%  de Energia (EMS)

%% Registro de variáveis existente no workspace antes da execução do script

% Registra variáveis existentes no workspace neste instante
if (exist('previousVars', 'var'))
    previousVars{:, end + 1} = who;
else
    previousVars{1} = who;
end

%% Parâmetros de controle

rectifier.control.pwm.f_s = 10e+3;  % Frequência de chaveamento dos PWMs de controle do circuito de retificador [Hz]
rectifier.control.pwm.v_clear = 0;  % Tensão de 'clear' dos PWMs de controle do circuito de retificador [V]
rectifier.control.pwm.v_set = 5;    % Tensão de 'set' dos PWMs de controle do circuito de retificador [V]
rectifier.control.pwm.v_out = 15;   % Tensão de saída do driver dos PWMs de controle do circuito de retificador [V]

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
