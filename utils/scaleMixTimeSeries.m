%% Recebe séries temporais e fatores de escala para retornar uma composição
 % delas escaladas

function scaledMixedTimeSeries = scaleMixTimeSeries(varargin)

% Inicialização da variável de retorno
scaledMixedTimeSeries = NaN;

% Verificação da quantidade de argumentos de entrada
if (mod(nargin, 2) ~= 0)
    disp('Número de parâmetros de entrada incorreto');
    return;
end

% Inicializações necessárias ao funcionamento do laço
timeSeriesMat = [];

% Laço de execução
for cntr = 1:(nargin/2)
    % Indexadores dos parâmetros de entrada
    indexTS = (cntr - 1)*2 + 1;         % Indexador das séries temporais
    indexSF = (cntr - 1)*2 + 2;         % Indexador dos fatores de escala
    
    % Variáveis temporárias contendo os parâmetros de interesse para essa
    % iteração
    timeSeries = varargin{indexTS};     % Série temporal
    scaleFactor = varargin{indexSF};    % Fator de escala
    
    % Conversão da série temporal para forma matricial já aplicando fator
    % de escala ao campo de dados
    tmpTimeSeriesMat = [];
    tmpTimeSeriesMat(:, 1) = timeSeries.Time;
    tmpTimeSeriesMat(:, 2) = timeSeries.Data*scaleFactor;
    
    % Concatenação da série temporal em forma matricial desta iteração com
    % aquelas de iterações passadas
    timeSeriesMat = [timeSeriesMat; tmpTimeSeriesMat];
end

% Ordenação da série temporal em forma matricial
timeSeriesMat = sortrows(timeSeriesMat, 1);

% Conversão da forma matricial para objeto série temporal
time = timeSeriesMat(:, 1);
data = timeSeriesMat(:, 2);
scaledMixedTimeSeries = timeseries(data, time);

end

