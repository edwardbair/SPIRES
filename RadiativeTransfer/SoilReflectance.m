function [ R ] = SoilReflectance( lambda, soilType, varargin )
% [ R ] = SoilReflectance( lambda, soilType [,units] )
% spectral reflectance of soils
%
% INPUT
%   lambda - wavelength(s), default in um
%   soilType - either 'DarkLoam', 'BrownLoam', 'GraySiltyLoam',
%       'PaleBrownSilty', or 'DarkBrownSilty'
%
% OPTIONAL INPUT
%   units for wavelength, default 'um', but 'nm', 'mum', 'mm', 'cm', 'm', 'GHz'
%       are supported
%
% OUTPUT
%   R - reflectance

persistent already soilF minmaxWave

possibleSoilType = {'DarkLoam', 'BrownLoam', 'GraySiltyLoam',...
    'PaleBrownSilty', 'DarkBrownSilty'};

warning('deprecated function, use ECOSTRESS_SpecRefl instead')

p = inputParser;
% inputs and defaults
defaultUnits = 'um';
nonnegativeFcn = @(x) isnumeric(x) && all(x>=0);
addRequired(p,'lambda',nonnegativeFcn);
addRequired(p,'soilType',@ischar);
addOptional(p,'units',defaultUnits,@ischar);
parse(p,lambda,soilType,varargin{:});
% 1st pass, calculate the interpolation functions
if isempty(already)
    already = true;
    load ('SoilReflectanceSolar.mat') %#ok<LOAD>
    wvl = convertUnits(soilT.wvl,soilT.Properties.VariableUnits{1},defaultUnits);
    minmaxWave = [min(wvl) max(wvl)];
    for k=1:length(possibleSoilType)
        F = fit(wvl,soilT.(possibleSoilType{k}),'pchipinterp');
        soilF{k} = F;
    end
end

% all passes, do the interpolation
% if lambda is a matrix, need to vectorize and then reshape results
waveSize = size(lambda);
lambda = lambda(:);
% identify soil Type in the function table
found = false;
for k=1:length(possibleSoilType)
    if strcmpi(soilType,possibleSoilType{k})
        F = soilF{k};
        found = true;
        break
    end
end

if found
    % interpolate, setting 0 outside of range
    lambda = convertUnits(lambda,p.Results.units,defaultUnits);
    R = F(lambda);
    t = lambda<minmaxWave(1) | lambda>minmaxWave(2);
    if nnz(t)
        R(t) = 0;
    end
    if ~isequal(size(R),waveSize)
        R = reshape(R,waveSize);
    end
elseif strcmpi(soilType,'help')
    str = possibleSoilType{1};
    for k=2:length(possibleSoilType)
        str = cat(2,str,', ',possibleSoilType{k});
    end
    error('available soilTypes are %s',str)
else
    error(['soilType ' '''' soilType ''''...
        ' not recognized, use ''help'' to see list'])
end

end