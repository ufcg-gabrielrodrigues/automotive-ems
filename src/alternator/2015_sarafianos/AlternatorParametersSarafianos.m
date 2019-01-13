%% Características relativos ao alternador (Bosch LiX 180 A, 2.4 kW) 
%  utilizado no artigo de Sarafianos (2015): "Characterisation and
%  Modelling of Automotive Lundell Alternators"

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
    %% Realiza nova análise dos dados explícitos e implícitos disponíveis
    %  em Sarafianos (2015) e disponibiliza características do alternador
    
    % Execução dos scripts de análise
    ExplicitAlternatorParametersSarafianos;
    ImplicitAlternatorParametersSarafianos;
    
    % Carrega arquivos .MAT contendo características do alternador
    load(explicitAlternatorCharacteristics);
    load(implicitAlternatorCharacteristics);
    
    %% Registro de parâmetros definidos em estrutura que representa o alternador
    
    alternator.frictionWindageLosses = frictionWindageLosses;                   % Perdas por atrito e enrolamento [W]
    alternator.ironLoss = extFunctionHandle(@(n_r, i_f) n_r.*ironLoss(i_f));    % Perdas no ferro [W]
    alternator.n = n;                                                           % Número de fases
    alternator.p = p;                                                           % Número de pares de polos por fase
    alternator.rotor.l.value = 0.2;                                             % Indutância própria do circuito de excitação [H]
    alternator.rotor.r.value = r_f;                                             % Resistência do circuito de excitação a 20oC [Ohm]
    alternator.rotor.s.v_gs_on = 10;                                            % Tensão gate-source para acionamento do MOSFET do circuito de excitação [V]
    alternator.rotor.control.pwm.f_s = 10e+3;                                   % Frequência de chaveamento do PWM de controle do circuito de excitação [Hz]
    alternator.stator.connection = connection;                                  % Tipo de conexão do circuito de estator
    alternator.stator.slots = slots;                                            % Número de ranhuras no estator
    alternator.k_e.function = extFunctionHandle(@(i_f) ...
        (openCircuitVoltage(i_f)./i_f).*(1./(p.*pi./30)));                      % Constante de acoplamento el�trico [V/((rad/s)*A)]
    alternator.stator.l.function = inductance;                                  % Indutância própria por fase do circuito de armadura [H]
    alternator.stator.r.value = r_a;                                            % Resistência por fase do circuito de armadura a 20oC [Ohm]
    
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
