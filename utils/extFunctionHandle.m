%% Converte manipuladores de função para sua versão expandida e retorna no
 % tipo desejado

function extFun = extFunctionHandle(fun, varargin)

% Converte manipulador de função para expressão simbólica. Nessa conversão
% a expressão é expandida.
extFun = sym(fun);

% Verifica a presença de parâmetros opcionais
if (nargin > 1)
    % Salva o parâmetro opcional em variável que identifica o tipo de
    % retorno desejado (manipulador de função, expressão simbólica ou
    % 'string')
    type = varargin{1};
    
    if (strcmp(type, 'handle'))
        extFun = matlabFunction(extFun);
    elseif (strcmp(type, 'sym'))
        % Já é o tipo atual da variável
    elseif (strcmp(type, 'string'))
        extFun = matlabFunction(extFun);
        extFun = func2str(extFun);
    else
        extFun = matlabFunction(extFun);
    end
else
    extFun = matlabFunction(extFun);
end

end

