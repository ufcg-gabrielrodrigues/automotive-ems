%% Versão multivariável da função 'strcmp' que realiza comparação de 'string'

function logicalOut = multStrCmp(strBase, strList)

% Dimensão da base de busca
[r, c] = size(strBase);

% Dimensão da lista de elementos buscados
n = length(strList);

% Inicialização da matriz lógica de saída
logicalOut = false(r, c);

% Laço de busca e atribuição
for i = 1:n
    % Busca pelo i-ésimo elemento da lista
    index = strcmp(strBase, strList{i});
    
    % Atualização do vetor lógico de saída
    logicalOut = logicalOut|index;
end

end

