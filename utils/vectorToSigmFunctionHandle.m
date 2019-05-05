%% Recebe vetor de parâmetros e retorna manipulador para função sigmoide

function sigmFunctionHandle = vectorToSigmFunctionHandle(v, varargin)

% Registra parâmetro opcional 'nome da variável'
if (nargin > 1)
    var = varargin{1};
else
    var = 'x';
end

% Determina quantidade de parâmetros e, consequentemente, ordem da função
n = length(v);

% Processo montagem da função
if (length(v) == 4)
    sigmFunctionHandle = strcat('@(', var, ') ', num2str(v(1)), ' + (', ...
        num2str(v(2)), ' - ', num2str(v(1)), ')./(1 + 10.^((', ...
        num2str(v(3)), ' - ', var, ').*', num2str(v(4)), '))');
else
    disp('The dimension of parameters vector must be equal to 4');
    sigmFunctionHandle = @(x) NaN;
end

% Conversão do conjunto de caracteres representando a função para o
% manipulador para a função de fato
sigmFunctionHandle = str2func(sigmFunctionHandle);

end

