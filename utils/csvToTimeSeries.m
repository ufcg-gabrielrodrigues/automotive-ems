%% Retorna objeto do tipo 'timeseries' criado a partir de arquivo .CSV

function timeSeries = csvToTimeSeries(filePath)

% Leitura do arquivo .CSV
tmp = csvread(filePath);

% Criacao do objeto a ser retornado
timeSeries = timeseries(tmp(:, 2), tmp(:, 1));

end
