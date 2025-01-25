%% Realiza ajuste de curva de série temporal

function fittedTimeSeries = fitTimeSeries(timeSeries, fitType)

if (strcmp(fitType, 'sigmoid'))
    fittedTimeSeries = sigm_fit(timeSeries.Time, timeSeries.Data, [], [], false);
else
    param = fit(timeSeries.Time, timeSeries.Data, fitType);
    fittedTimeSeries = coeffvalues(param);
end

end
