%% Recebe vetor de parâmetros e retorna manipulador para função polinomial

function polyFunctionHandle = vectorToPolyFunctionHandle(v, varargin)

% Registra parâmetro opcional 'nome da variável'
if (nargin > 1)
    var = varargin{1};
else
    var = 'x';
end

% Determina quantidade de parâmetros e, consequentemente, ordem da função
n = length(v);

% Processo iterativo de montagem da função
polyFunctionHandle = strcat('@(', var, ') ');

for cntr = 1:n
    polyFunctionHandle = strcat(polyFunctionHandle, num2str(v(cntr)), '.*', ...
        var, '.^', num2str(n - cntr));
    
    if (cntr ~= n)
        polyFunctionHandle = strcat(polyFunctionHandle, ' + ');
    end
end

% Conversão do conjunto de caracteres representando a função para o
% manipulador para a função de fato
polyFunctionHandle = str2func(polyFunctionHandle);

end

