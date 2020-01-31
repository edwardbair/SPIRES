function outS = ScaledSnowAbsorptionModel(cosZ,effectiveRadius,varargin)
% outS = ScaledSnowAbsorptionModel(cosZ,effectiveRadius,...)
%ScaledSnowAbsorptionModel -- calculates parameters of spectral snow
%absorption features, around 1030 and 1240 nm for grain size, up to 800 nm
%for dust or soot
%
%Input
%   cosZ - cosine of illumination angle
%   effectiveRadius - default units are um

%Optional input, name-value pairs
%   radiusUnits - default 'um', but 'mm' and 'nm' and other units are supported
%   contam - 'dust' or 'soot'
%   contamRadius - default from SnowCloudLimits.m sqrt of max, same units
%       as effectiveRadius
%   contamConc - particulate concentration, dimensionless, default middle
%       of range from SnowCloudLimits.m
%   lookup - use Mie lookup tables, default true
%
%Output structure, 3 substructures corresponding to dust/soot, nm1030, and
%nm1240, each containing
%   name - name of feature
%   wavelength across feature
%   reflDiff - difference continuum minus actual, scaled to continuum
%       values, corresponding to wavelengths
%   continuumR - continuum reflectance, corresponding to wavelengths
%   actualR - actual reflectance, corresponding to wavelengths

%
%based on Nolin & Dozier, RSE 2000, DOI: 10.1016/S0034-4257(00)00111-5

p = inputParser;
rangeFcn = @(x) isnumeric(x) && all(x(:)<=1) && all(x(:)>0);
positiveFcn = @(x) isnumeric(x) && all(x(:)>0);
addRequired(p,'cosZ',rangeFcn)
addRequired(p,'effectiveRadius',positiveFcn)
addParameter(p,'radiusUnits','um',@ischar)
addParameter(p,'contam','neither',@ischar)
addParameter(p,'contamradius',[],positiveFcn)
addParameter(p,'contamconc',[],positiveFcn)
addParameter(p,'lookup',true,@islogical)
parse(p,cosZ,effectiveRadius,varargin{:})

% wavelengths that cover the absorption bands
allWave = (400:5:1370)';

% reflectance across the band
cleanSnow = strcmpi(p.Results.contam,'neither');
if cleanSnow
    R = SnowCloudSpectralRefl('snow','cosZ',p.Results.cosZ,'radius',p.Results.effectiveRadius,...
        'sizeUnit',p.Results.radiusUnits,'wavelength',allWave,'waveu','nm','lookup',p.Results.lookup);
    
    outS = ScaledSnowAbsorptionActual(allWave,'nm',R.refl,'lookup',p.Results.lookup);

else
    R = SnowCloudSpectralReflectance(p.Results.cosZ,p.Results.effectiveRadius,...
        p.Results.radiusUnits,allWave,'nm','contam',p.Results.contam,'contamconc',...
        p.Results.contamconc,'contamradius',p.Results.contamradius,...
        'lookup',p.Results.lookup);
    
    outS = ScaledSnowAbsorptionActual(allWave,'nm',R,...
        'snowRadius',p.Results.effectiveRadius,...
        'cosZ',p.Results.cosZ,'radiusUnits',p.Results.radiusUnits,...
        'lookup',p.Results.lookup);
    
end
end