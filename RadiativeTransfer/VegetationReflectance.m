function [ R ] = VegetationReflectance( lambda, vegType, varargin )
% [ R ] = VegetationReflectance( lambda, vegType [,units] )
% spectral reflectance of vegetation
% Data for the solar spectrum are from Rob Green.
% Data for the thermal IR spectrum are from Zhengming Wan's Emissivity
% Library https://icess.eri.ucsb.edu/modis/EMIS/html/em.html
%
% INPUT
%   wave - wavelength(s), default in um
%   vegType - either 'Conifer', 'Grass', 'BroadLeaf', 'SageBrush', or 'NPV'
%
% OPTIONAL INPUT
%   units for wavelength, default 'um', but 'nm', 'mum', 'mm', 'cm', 'm', 'GHz'
%       are supported
%
% OUTPUT
%   R - reflectance

persistent already vegF minmaxWave

possibleVegetationType = {'Conifer','Grass','BroadLeaf','SageBrush','NPV'};

warning('deprecated function, use ECOSTRESS_SpecRefl instead')

p = inputParser;
% inputs and defaults
defaultUnits = 'um';
nonnegativeFcn = @(x) isnumeric(x) && all(x>=0);
addRequired(p,'lambda',nonnegativeFcn);
addRequired(p,'vegType',@ischar);
addOptional(p,'units',defaultUnits,@ischar);
parse(p,lambda,vegType,varargin{:});
% 1st pass, calculate the interpolation functions
if isempty(already)
    already = true;
    load ('VegetationReflectanceSolar.mat') %#ok<LOAD>
    load ('VegetationReflectanceIR.mat') %#ok<LOAD>
    % combine the tables
    vegT = [vegT; WanVegT]; %#ok<NODEF>
    wvl = convertUnits(vegT.wvl,vegT.Properties.VariableUnits{1},defaultUnits);
    minmaxWave = [min(wvl) max(wvl)];
    for k=1:length(possibleVegetationType)
        F = fit(wvl,vegT.(possibleVegetationType{k}),'pchipinterp');
        vegF{k} = F;
    end
end

% all passes, do the interpolation
% if lambda is a matrix, need to vectorize and then reshape results
waveSize = size(lambda);
lambda = lambda(:);
% identify substance in the function table
found = false;
for k=1:length(possibleVegetationType)
    if strcmpi(vegType,possibleVegetationType{k})
        F = vegF{k};
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
elseif strcmpi(vegType,'help')
    str = possibleVegetationType{1};
    for k=2:length(possibleVegetationType)
        str = cat(2,str,', ',possibleVegetationType{k});
    end
    error('available vegTypes are %s',str)
else
    error(['vegType ' '''' vegType ''''...
        ' not recognized, use ''help'' to see list'])
end

end