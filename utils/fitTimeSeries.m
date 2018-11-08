%% Realiza ajuste de curva de série temporal

function fittedTimeSeries = fitTimeSeries(timeSeries, fitType)

fittedTimeSeries = fit(timeSeries.Time, timeSeries.Data, fitType);

end
