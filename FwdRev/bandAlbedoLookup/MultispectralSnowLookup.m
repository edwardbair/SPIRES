function [bandRefl,varargout] = MultispectralSnowLookup(sensor,varargin)
% [bandRefl] = MultispectralSnowLookup(sensor,Pscript and/or Name-Value pairs)
% [bandRefl,Pscript] = MultispectralSnowLookup(sensor,Pscript and/or Name-Value pairs)
% [bandRefl,Pscript,Tbl] = MultispectralSnowLookup(sensor,Pscript and/or Name-Value pairs)
%
%Calculates the band-integrated reflectance for several multispectral sensors,
%based on snow physical properties and illumination conditions
%
%Input data may consist of a prescription structure with the following
%fields, or the input can consist of Name/Value pairs, or a prescription
%can be modified by including the prescription and Name/Value pairs that
%modify it. If a prescription is entered, it must be the second argument
%(the first optional argument).
%
%Variables entered can be scalars, vectors, matrices, or N-dimensional
%objects, but any non-scalar variables must have the same size, and matrices
%or N-dimensional objects are converted to column vectors. To consider all
%possible combinations, use ndgrid.
%
%Required input
%   sensor - e.g. 'Landsat' or 'MODIS', case-insensitive
%Optional input, either as fields in the prescription in any order or as
%Name/Value pairs in any order are:
%Variables needed for clean or dirty snow
%   elevation - in meters, default 3000
%   mu0 - cosine of the solar zenith angle on a horizontal surface
%   muS - cosine of the illumination angle on a sloping surface
%   radius or SSA - radius of the snow grains in um, or specific surface
%       area in m2/kg (either can be entered, but not both)
%Variables needed only for dirty snow
%   dust - concentration as a mass fraction
%   soot - concentration as a mass fraction
%
%Output
%   bandRefl - reflectance of the multispectral bands as a NxBands matrix,
%       where N is the size of the input variables (N=1 if all scalar) and
%       Bands is the # of bands
%Optional output
%   Pscript - prescription as used in the calculation, can be modified for
%       different values of the variables
%   Tbl - which includes the input values as well as the bandRefl values

persistent FClean FDirty SClean SDirty lastSensor

S = SnowCloudLimits;
defaultElevation = 3000;
p = inputParser;
addRequired(p,'sensor',@ischar)
if isstruct(varargin{1})
    addOptional(p,'Pscript',struct([]),@isstruct)
end
addParameter(p,'elevation',3000,@isnumeric)
addParameter(p,'mu0',[],@(x) isnumeric(x) && (all(x(:)>0) && all(x(:)<=1)))
addParameter(p,'muS',[],@(x) isnumeric(x) && (all(x(:)>=0) && all(x(:)<=1)))
addParameter(p,'radius',[],@(x) isnumeric(x) &&...
    (all(x(:)>=.99*S.snowRadius(1)) && all(x(:)<=1.01*S.snowRadius(2))))
addParameter(p,'ssa',[],@(x) isnumeric(x) &&...
    (all(x(:)>=.99*S.snowSSA(1)) && all(x(:)<=1.01*S.snowSSA(2))))
addParameter(p,'dust',0,@(x) isnumeric(x) &&...
    (all(x(:)>=0) && all(x(:)<=1.01*S.dust(2))))
addParameter(p,'soot',0,@(x) isnumeric(x) &&...
    (all(x(:)>=0) && all(x(:)<=1.01*S.soot(2))))
parse(p,sensor,varargin{:})

%% parse the inputs
if isfield(p.Results,'Pscript')
    Pscript = p.Results.Pscript;
else
    Pscript = struct([]);
end
possibleFields = {'elevation','mu0','muS','ssa','radius','dust','soot'};
assert(isempty(p.Results.radius) || isempty(p.Results.ssa),...
    'either ''radius'' or ''SSA'' can be entered, but not both')
if ~isempty(Pscript)
    % replace any that are in the Name/Value pairs
    for k=1:length(possibleFields)
        if ~isempty(p.Results.(possibleFields{k}))
            if strcmpi(possibleFields{k},'ssa')
                Pscript.SSA = p.Results.ssa;
            elseif strcmpi(possibleFields{k},'radius')
                Pscript.SSA = radius2SSA(p.Results.radius,'um');
            else
                Pscript.(possibleFields{k}) = p.Results.(possibleFields{k});
            end
        end
    end
else
    Pscript = struct;
    % add any that are in the Name/Value pairs
    for k=1:length(possibleFields)
        if ~isempty(p.Results.(possibleFields{k}))
            if strcmpi(possibleFields{k},'ssa')
                Pscript.SSA = p.Results.ssa;
            elseif strcmpi(possibleFields{k},'radius')
                Pscript.SSA = radius2SSA(p.Results.radius,'um');
            else
                Pscript.(possibleFields{k}) = p.Results.(possibleFields{k});
            end
        end
    end
end
% check to make sure we have the 4 essential variables to calculate
% reflectance of the bands
if ~isfield(Pscript,'elevation')
    Pscript.elevation = defaultElevation;
end
assert(isfield(Pscript,'SSA'),'either ''radius'' or ''SSA'' must be provided')
assert(isfield(Pscript,'mu0'),'''mu0'' must be provided')
if ~isfield(Pscript,'muS')
    Pscript.muS = Pscript.mu0;
end
%% sizes of the inputs
[elevation,mu0,muS,SSA,dust,soot] = checkSizes(Pscript.elevation,...
    Pscript.mu0,Pscript.muS,Pscript.SSA,Pscript.dust,Pscript.soot);
elevation = elevation(:);
mu0 = mu0(:);
muS = muS(:);
SSA = SSA(:);
dust = dust(:);
soot = soot(:);
%% identify the sensor
assert(~isempty(sensor),'first argument, sensor name, must be specified')
if contains(sensor,'landsat','IgnoreCase',true)
    sensor = 'LandsatOLI';
elseif contains(sensor,'modis','IgnoreCase',true)
    sensor = 'MODIS';
else % others to be added
    error('sensor ''%s'' not recognized',sensor)
end

%% load the lookup functions
if isempty(FClean) || isempty(lastSensor) || ~strcmp(lastSensor,sensor)
    functionFile = [sensor '_Lookup.mat'];
    load(functionFile,'FClean','SClean','FDirty','SDirty')
    assert(isequal(SDirty.sensor,SClean.sensor),...
        'clean and dirty lookup functions are not for the same sensor')
    assert(isequal(SDirty.bands,SClean.bands),...
        'clean and dirty lookup functions do not include the same bands')
end

bandRefl = zeros(length(elevation),length(SClean.bands));
doDirt = any(soot>0) || any(dust>0);
for b=1:length(SClean.bands)
    scaleSSA = bandScaling(sensor,SClean.bands(b),'SSA',SSA);
    cleanRefl = FClean{b}(elevation,mu0,scaleSSA,muS);
    if doDirt
        scaleDust = bandScaling(sensor,SDirty.bands(b),'dust',dust);
        scaleSoot = bandScaling(sensor,SDirty.bands(b),'soot',soot);
        dirtyDiff = FDirty{b}(elevation,mu0,scaleSSA,muS,...
            scaleDust,scaleSoot);
        bandRefl(:,b) = cleanRefl+dirtyDiff;
    else
        bandRefl(:,b) = cleanRefl;
    end
end

if nargout>1
    fields = fieldnames(Pscript);
    for k=1:length(fields)
        Pscript.(fields{k}) = unique(Pscript.(fields{k}));
    end
    varargout{1} = Pscript;
    if nargout>2
        varargout{2} = table(elevation,mu0,SSA,muS,dust,soot,bandRefl);
    end
end
lastSensor = sensor;
end