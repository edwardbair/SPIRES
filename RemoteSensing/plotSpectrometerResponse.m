function plotSpectrometerResponse(varargin)
%h = plotSpectrometerResponse([sensor - default AVIRIS-NG])
%
p = inputParser;
addOptional(p,'sensor','AVIRIS-NG',@ischar)
parse(p,varargin{:})

units = 'nm';
T = SensorTable(p.Results.sensor,units);
% calculate and plot, based on FWHM
for k=1:height(T)
    meanValue = T.CentralWavelength(k);
    lowerValue = T.LowerWavelength(k);
    sighat = fzero(@fwhm,T.Bandwidth(k)/2);
    wave = linspace(T.LowerWavelength(k)-T.Bandwidth(k)/10,...
        T.UpperWavelength(k)+T.Bandwidth(k)/10,50);
    thisPD = makedist('normal','mu',meanValue,'sigma',sighat);
    y = pdf(thisPD,wave);
    y = y/max(y);
    plot(wave,y,'k')
    hold on;
end

    function half = fwhm(sigma)
        pd = makedist('normal','mu',meanValue,'sigma',sigma);
        half = pdf(pd,lowerValue)/pdf(pd,meanValue)-0.5;
    end
end