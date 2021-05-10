function [ R ] = SnowCloudSensorReflectance(cosZ,radius,radiusUnits,sensor,bands, varargin)
% [ R ] = SnowCloudSensorReflectance(cosZ,radius,radiusUnits,sensor,bands [, 'R0', R0, 'WE', value, 'weUnits', '(units)', 'substance', '(ice or water)'])
%snow or cloud reflectance integrated over sensor bands
%
%Inputs
% cosZ - cosine of illumination angle (same size as radius)
% radius, scalar or vector (same size as cosZ)
% radiusUnits - any metric length ('angstrom', 'nm', 'um' or 'mum', 'mm', 'cm', 'm')
% sensor - character string identifying sensor in SensorTable.m
% bands - cell vector of character strings identifying bands or just
%   numbers if bands for this sensor are by number
%Optional input, name-value pair
% 'WE' - snow or cloud water equivalent (Inf if not specified)
% 'weUnits' - any metric length (typically 'mm', 'cm', or 'm')
% 'R0' - reflectance of underlying soil, either scalar or size of bands
%   vector (ignored unless WE specified)
% 'substance' - 'ice' (default) or 'water'
% 'contam', followed by 'soot' or 'dust' (program deals with either, but
%       not both in combination, default is 'neither', in which case snow is
%       assumed clean)
% 'contamSize', followed by effective radius of contaminant, same 'units'
%       as for grain size (default 2.5 um for dust, 300 nm for soot)
% 'contamConc', concentration of contaminants, as a mass fraction
% 'lookup', if true, lookup tables are called if they exist, with
%   a warning if they don't or if input values are out of range, if false,
%   calculations are done, with no warning
% 'ignoreSolar', if true ignores solar radiation and just provides
%   band-average reflectivity, default false, in which case solar radiation
%   accounted for unless outside range of SolarScale.m
%   (this is needed to calculate emissivities around 4 um)
%
%Output
% R - snow reflectance, dimensionless, size # radii x # bands
%

warning('function %s deprecated, use the toolbox/SnowCloudReflectance instead', mfilename)

% parse inputs
defaultWE = Inf;
defaultWEunits = '';
defaultR0 = 0;
defaultSubstance = 'ice';
defaultContam = 'neither';
validationFcn = @(x) isnumeric(x) || iscell(x);
rangeFcn = @(x) isnumeric(x) && all(x(:)>=0 & x(:)<=1);
positiveFcn = @(x) isnumeric(x) && all(x(:)>0);
nonnegativeFcn = @(x) isnumeric(x) && all(x(:)>=0);
p = inputParser;
addRequired(p,'cosZ',rangeFcn)
addRequired(p,'radius',@isnumeric)
addRequired(p,'radiusUnits',@ischar)
addRequired(p,'sensor',@ischar)
addRequired(p,'bands',validationFcn)
addParameter(p,'we',defaultWE,positiveFcn)
addParameter(p,'weunits',defaultWEunits,@ischar)
addParameter(p,'r0',defaultR0,rangeFcn)
addParameter(p,'substance',defaultSubstance,@ischar)
addParameter(p,'contam',defaultContam,@ischar)
addParameter(p,'contamsize',0,positiveFcn)
addParameter(p,'contamconc',[],nonnegativeFcn)
addParameter(p,'lookup',false,@islogical)
addParameter(p,'ignoresolar',false,@islogical)
parse(p,cosZ,radius,radiusUnits,sensor,bands,varargin{:})

% check sizes if not scalars
if strcmp(p.Results.contam,'neither')
    [cosZ,radius,WE] = checkSizes(p.Results.cosZ,p.Results.radius,p.Results.we);
    contamConc = zeros(size(cosZ));
else
    [cosZ,radius,WE,contamConc] =...
        checkSizes(p.Results.cosZ,p.Results.radius,p.Results.we,p.Results.contamconc);
end
R0 = p.Results.r0;
% convert band numbers to cell vectors if input is numeric
if isnumeric(p.Results.bands)
    band = cell(size(p.Results.bands));
    for k=1:length(band)
        band{k} = num2str(p.Results.bands(k));
    end
else
    band = p.Results.bands;
end

if isscalar(R0) && ~isscalar(band)
    R0 = repmat(R0,size(band));
elseif ~isscalar(band) && ~isscalar(R0)
    assert(isequal(size(R0),size(band)),...
        'if not scalars, sizes of R0 and bands must be same')
end

if p.Results.lookup
    %copy revised input into input structure
    S = p.Results;
    S.cosZ = cosZ;
    S.radius = radius;
    S.bands = band;
    S.contamconc = contamConc;
    S.r0 = R0;
    S.we = WE;
    S.sensor = lower(p.Results.sensor);
    [R,flag] = lookupSensor(S);
    if isempty(flag)
        return
    else
        warning('lookup table failed, flag=''%s'', continuing with Mie calculations',flag)
    end
end

% table for this sensor and its bands
X = SensorTable(p.Results.sensor);
t = false(length(X.Band),1);
for k=1:length(band)
    t = t | X.Band==band{k};
end
X = X(t,:);

R = zeros(length(radius),height(X));
for k=1:height(X)
    for n=1:length(radius)
        R(n,k) = SnowCloudBandPassReflectance(cosZ(n),radius(n),p.Results.radiusUnits,...
            [X{k,'LowerWavelength'} X{k,'UpperWavelength'}],'um','R0',R0(k),...
            'WE',WE(n),'weUnits',p.Results.weunits,'substance',...
            p.Results.substance,'contam',p.Results.contam,...
            'contamSize',p.Results.contamsize,'contamConc',contamConc(n),...
            'ignoreSolar',p.Results.ignoresolar);
    end
end

end

function [R,flag] = lookupSensor(Results)
persistent PreviousSensorForTable LUT
flag = '';
switch Results.substance
    case 'water'
        SensorForTable = ['LUT_' Results.sensor '_WaterCloud'];
        nVar = 4;
    case 'ice'
        switch min(Results.we)
            case Inf % snow
                nVar = 3;
                switch Results.contam
                    case 'neither' %use dusty snow table with contamConc=0
                        SensorForTable = ['LUT_' Results.sensor '_SnowDust'];
                    case 'dust'
                        SensorForTable = ['LUT_' Results.sensor '_SnowDust'];
                    case 'soot'
                        SensorForTable = ['LUT_' Results.sensor '_SnowSoot'];
                    otherwise
                        error('''contam'' ''%s'' not recognized',Results.contam)
                end
            otherwise % ice cloud
                nVar = 4;
                SensorForTable = ['LUT_' Results.sensor '_IceCloud'];
        end
    otherwise
        error('''substance'' ''%s'' not recognized',Results.substance)
end
if isempty(PreviousSensorForTable) ||...
        ~strcmp(PreviousSensorForTable,SensorForTable)
    if exist([SensorForTable '.mat'],'file')==2
        LUT = load([SensorForTable '.mat']);
        PreviousSensorForTable = SensorForTable;
    else
        flag = ['file ' SensorForTable '.mat' ' not found'];
        R = [];
        return
    end
end
%make sure for the right sensor
if ~strcmp(LUT.sensor,Results.sensor)
    flag = [SensorForTable ' is for sensor ' LUT.sensor, ' not ' Results.sensor];
    R = [];
    return
end

% units
Results.radius = convertUnits(Results.radius,Results.radiusUnits,'mum');
if nVar==4
    Results.we = convertUnits(Results.we,Results.weunits,'mm');
end

R = zeros(length(Results.radius),length(Results.bands));
for k=1:length(Results.bands)
    t = strcmp(LUT.bands,Results.bands{k});
    n = find(t);
    assert(length(n)==1,'%d matches to LUT.bands for band{%d}=%s',...
        length(n),k,Results.bands{k})
    F = LUT.F{n}; %gridded interpolant
    switch nVar
        case 3
            R(:,k) = F(Results.cosZ,Results.radius,Results.contamconc);
        case 4
            R(:,k) = F(Results.cosZ,Results.radius,Results.we,Results.r0(k));
    end
    if any(isnan(R(:,k)))
        f = find(isnan(R(:,k)),1,'first');
        switch nVar
            case 3
                flag = sprintf('%s: some input values out of range, first bad value cosZ=%g, radius=%g, contamConc=%g',...
                    SensorForTable,Results.cosZ(f),Results.radius(f),Results.contamconc(f));
            case 4
                flag = sprintf('%s: some input values out of range, first bad value cosZ=%g, radius=%g, WE=%g, R0=%g',...
                    SensorForTable,Results.cosZ(f),Results.radius(f),Results.we(f),Results.r0(f));
        end
        R = [];
        return
    end
end
end