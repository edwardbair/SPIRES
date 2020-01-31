function [h] = plotSensorSpec(sensor,varargin)
%plot specs for a sensor in SensorTable
%
%Input
%   sensor - character vector of sensor
%Optional name-value pairs
%   'bands' - numeric or categorical vector, default all
%   'units' - typically either 'um' or 'nm' but any length unit works,
%       default 'um'
%   'scaleX' - default 'log', any alternative makes linear
%   'scaleY' - default 'log', any alternative makes linear
%   'color' - either letter or 3-element vector
%   'linewidth' - width of lines (pts), scalar

p = inputParser;
addRequired(p,'sensor',@ischar)
addParameter(p,'bands',[],@(x) isnumeric(x) || iscategorical(x))
addParameter(p,'units','um',@ischar)
addParameter(p,'scaleX','log',@ischar)
addParameter(p,'scaleY','log',@ischar)
addParameter(p,'color','',@(x) isnumeric(x) || ischar(x))
addParameter(p,'linewidth',20,@isnumeric)

parse(p,sensor,varargin{:})

%input values, will fail if sensor not supported
T = SensorTable(sensor,p.Results.units);

% plot specifications
if isempty(p.Results.bands)
    for k=1:height(T)
        h = plotrow(T(k,:));
    end
else
    if iscategorical(p.Results.bands)
        bands = p.Results.bands;
    else
        bands = categorical(p.Results.bands);
    end
    for b=1:length(bands)
        t = T.Band==bands(b);
        h = plotrow(T(t,:));
    end
end
set(gca,'YDir','reverse')

    function h = plotrow(thisRow)
        if strcmpi(p.Results.scaleX,'log') && strcmpi(p.Results.scaleY,'log')
            plotFun = @loglog;
        elseif strcmpi(p.Results.scaleX,'log')
            plotFun = @semilogx;
        elseif strcmpi(p.Results.scaleY,'log')
            plotFun = @semilogy;
        else
            plotFun = @plot;
        end
        x = [thisRow.LowerWavelength thisRow.UpperWavelength];
        y = [thisRow.SpatialResolution thisRow.SpatialResolution];
        if isempty(p.Results.color)
            h = plotFun(x,y,'LineWidth',p.Results.linewidth);
        else
            h = plotFun(x,y,'LineWidth',p.Results.linewidth,'Color',p.Results.color);
        end
        hold on;
    end
end

