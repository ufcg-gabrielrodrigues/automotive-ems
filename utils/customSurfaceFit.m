function [fitresult, gof] = customSurfaceFit(x, y, z, f, fx, fy, fz, p0, varargin)
    % 
    if (nargin > 8)
        plot_flag = varargin{1};
    else
        plot_flag = false;
    end
    
    % 
    [xData, yData, zData] = prepareSurfaceData(x, y, z);

    % 
    ft = fittype(f, 'independent', {fx, fy}, 'dependent', fz);
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.StartPoint = p0;

    % 
    [fitresult, gof] = fit([xData, yData], zData, ft, opts);

    if (plot_flag == true)
        % 
        figure('Name', 'Custom Surface Fit');
        h = plot( fitresult, [xData, yData], zData );
        legend( h, 'Custom Surface Fit', 'z vs. x, y', 'Location', 'NorthEast' );
        grid on
    end
end