function [ Tb ] = SnowCloudSensorBrightnessTemp(physTemp,cosZ,sensor, bands, varargin)
% [ Tb ] = SnowCloudSensorBrightnessTemp(physTemp,cosZ,sensor, bands, varargin)
%brightness temperature of snow or cloud integrated over wavelength range
%
%Inputs
% physTemp - physical temperature, Kelvin, scalar
% cosZ - cosine of viewing angle, scalar
% sensor - character string identifying sensor in SensorTable.m
% bands - cell vector of character strings or numeric vector of band numbers
%Optional input, name-value pair
% 'radius' scalar
% 'radiusUnits' - any metric length ('angstrom', 'nm', 'um' (default) or 'mum', 'mm', 'cm', 'm')
% 'WE' - snow or cloud water equivalent (Inf if not specified)
% 'weUnits' - any metric length (typically 'mm', 'cm', or 'm')
% 'substance' - 'ice' (default) or 'water'
%
%Output
% Tb - brightness temperature, Kelvin, same size as bands vector

% parse inputs
defaultWE = Inf;
defaultWEunits = '';
defaultSubstance = 'ice';
defaultRadius = 250;
defaultRadiusUnits = 'um';
p = inputParser;
validationFcn = @(x) isnumeric(x) || iscell(x);
rangeFcn = @(x) isnumeric(x) && (all(x>=0 & x<=1));
positiveFcn = @(x) isnumeric(x) && (all(x>0));
addRequired(p,'physTemp',positiveFcn)
addRequired(p,'cosZ',rangeFcn)
addRequired(p,'sensor',@ischar)
addRequired(p,'bands',validationFcn)
addParameter(p,'radius',defaultRadius,positiveFcn)
addParameter(p,'radiusunits',defaultRadiusUnits,@ischar)
addParameter(p,'we',defaultWE,positiveFcn)
addParameter(p,'weunits',defaultWEunits,@ischar)
addParameter(p,'substance',defaultSubstance,@ischar)
parse(p,physTemp,cosZ,sensor,bands,varargin{:})

% table for this sensor and its bands
assert(ismatrix(p.Results.bands) &&...
    (size(p.Results.bands,1)==1 || size(p.Results.bands,2)==1),...
    'bands must be a vector (row or column), not a 2D or multi-D object')
% band a cell vector of strings (some sensors use alphanumeric
% designators)
if isnumeric(p.Results.bands)
    band = cell(size(p.Results.bands));
    for k=1:length(p.Results.bands)
        band{k} = num2str(p.Results.bands(k));
    end
else
    band = p.Results.bands;
end
X = SensorTable(p.Results.sensor);
t = false(length(X.Band),1);
for k=1:length(band)
    t = t | strcmpi(X.Band,band{k});
end
X = X(t,:);

Tb = zeros(size(band));
for k=1:height(X)
    Tb(k) = SnowCloudBandPassBrightnessTemp(physTemp,cosZ,...
        [X{k,'LowerWavelength'},X{k,'UpperWavelength'}],'mum',...
        'radius',p.Results.radius,'radiusUnits',p.Results.radiusunits,...
        'WE',p.Results.we,'weUnits',p.Results.weunits,...
        'substance',p.Results.substance);
end
end