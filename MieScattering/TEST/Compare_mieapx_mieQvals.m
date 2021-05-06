function [Tbl] = Compare_mieapx_mieQvals(substance,wavelength,radius)
%compare calculations from mieapx and mieQvals
%
%input
% substance - 'ice' or 'water'
% wavelength and radius must have same size if not scalars and must be in
%   the same units (um)
%
%output
% Tbl of comparisons, also presented as a plot

waveScalar = isscalar(wavelength);
sizeScalar = isscalar(radius);
[wavelength,radius] = checkSizes(wavelength,radius);
CRefIn = RefractiveIndex(wavelength,substance);

xx = 2*pi*radius./wavelength;
Mapx = mieapx(CRefIn,xx);
MQ = mieQvals(CRefIn,xx);

Tbl = table(wavelength(:),radius(:),CRefIn(:),Mapx.omega(:),Mapx.g(:),...
    MQ.omega(:),MQ.g(:),...
    'VariableNames',{'wavelength','radius','RefIndex','omega_apx','g_apx',...
    'omega_Q','g_Q'});
subplot(1,2,1)
if sizeScalar
    loglog(wavelength,[1-Mapx.omega 1-MQ.omega],'linewidth',1);
    xlabel('log wavelength')
    ylabel('1-omega')
elseif waveScalar
    loglog(radius,[1-Mapx.omega 1-MQ.omega],'linewidth',1);
    xlabel('log radius')
    ylabel('1-omega')
else
    imagesc(unique(log10(radius(:))),unique(log10(wavelength(:))),...
        log10(1-(Mapx.omega-MQ.omega)))
    xlabel('log radius')
    ylabel('log wavelength')
    title('log10(1-(omega diff))')
end
subplot(1,2,2)
if sizeScalar
    semilogx(wavelength,[Mapx.g MQ.g],'linewidth',1)
    xlabel('log wavelength')
    ylabel('g')
elseif waveScalar
    semilogx(radius,[Mapx.g MQ.g],'linewidth',1)
    xlabel('log radius')
    ylabel('g')
else
    imagesc(unique(log10(radius(:))),unique(log10(wavelength(:))),...
        Mapx.g-MQ.g)
    xlabel('log radius')
    ylabel('log wavelength')
    title('g diff')
end

end