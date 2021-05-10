function [ R ] = SnowCloudReflectance(cosZ,substance,radius, radiusUnit, varargin)
% [ R ] = SnowCloudReflectance(cosZ,substance,radius, radiusUnit, varargin)
%reflectance of snow or cloud across band passes or sensor bands
%
%Inputs (if not scalars, must be of same size)
% cosZ - cosine of illumination angle
% substance - 'ice' or 'snow' (same), 'water', 'wetSnow', or
% 'mixedPhaseCloud'
%   (to model a mixed phase cloud, specifiy 'mixedPhaseCloud' as the
%   substance, iceRadius as the radius, 'waterConc' as the mass fraction of
%   water, and 'cloudRadius' as the size of the water droplets)
% radius - either scalar or vector
% radiusUnit - any metric length ('angstrom', 'nm', 'um' or 'mum', 'mm', 'cm', 'm')
%
%Optional input, name-value pairs in any order
%Arguments about the bands (either 'bandPass' or 'sensor'/'bands' must
%be specified, but not both
% 'bandPass' - matrix of size Nx2, specifying the wavelength ranges of N band
%   passes
% 'waveUnit' - wavelength units of bandPass (not needed when 'sensor' is
%   specified)
% 'sensor' - followed by character string, any sensor from SensorTable.m
% 'band' - either numeric vector, cell vector, or categorical vector of bands,
%   of, if omitted, all bands for that sensor
%   (numeric vector can be used for those sensors whose bands are designated just by number)
% 'ignoreSolar', if true ignores solar radiation and just provides
%   band-average reflectivity, default false, in which case solar radiation
%   accounted for unless outside range of SolarScale.m
%   (this is needed to calculate emissivities around 4 um)
%Arguments about the snow or cloud, either scalar or size of radius unless
%noted
% 'WE' - snow or cloud water equivalent (Inf if not specified)
% 'weUnit' - any metric length (typically 'mm', 'cm', or 'm'), default 'mm'
% 'waterConc' - mass fraction of water, needed if 'wetSnow' or
%   'mixedPhaseCloud' is the substance
% 'R0' - reflectance of underlying ground (ignored unless WE specified),
%   can be scalar or same size as number of bands, default 0.0,
% 'dust' - dust concentration (mass fraction)
% 'dustRadius' = optical radius of dust
% 'soot' - soot concentration (mass fraction)
% 'sootRadius' - optical radius of soot
% 'cloudRadius' - needed if substance is 'mixedPhaseCloud', same units as
%   radius
% 'lookup', if true, lookup tables are called for Mie calculations if they
%   exist, with a warning if they don't or if input values are out of range,
%   if false, Mie calculations are done (default true)
%
%Output
% R - snow or cloud reflectance, dimensionless, same size dimensions as
%   radius (or cosZ) x number of bands

% hang onto wavelength values corresponding to refractive indices
persistent wvice

warning('function %s deprecated, use the toolbox/SnowCloudReflectance instead', mfilename)

narginchk(5,25)
nargoutchk(0,1)

% parse inputs (convert radii to mum, wavelengths to nm, WE to mm,
optargin=size(varargin,2);
assert (mod(optargin,2)==0,'must be even number of optional arguments')
iStruct = parseInput(cosZ,substance,radius,radiusUnit,varargin{:});

% if spectrometer, call SnowCloudSpectralReflectance directly
if iStruct.spectrometer
    argc = varargin;
    nParam = length(varargin)/2;
    for k=nParam*2-1:-2:1
        if strcmpi(argc{k},'sensor') || strcmpi(argc{k},'bandPass')
            argc(k:k+1) = [];
        end
        if strcmpi(argc{k},'waveUnit')
            argc{k} = 'lambdaUnits';
        end
    end
    if isempty(argc)
        R = SnowCloudSpectralReflectance(cosZ,substance,radius,radiusUnit,...
            iStruct.wavelength,iStruct.waveUnit);
    else
        R = SnowCloudSpectralReflectance(cosZ,substance,radius,radiusUnit,...
            iStruct.wavelength,iStruct.waveUnit,argc{:});
    end
    return
end

% wavelengths to consider is the full range of all band passes at
% resolution of refractive index data
if isempty(wvice)
    [~,wvice] = RefractiveIndex([],'ice',iStruct.waveUnit);
end

%wavelengths to cover all band passes
k1 = find(wvice<min(iStruct.bandPass(:)),1,'last');
k2 = find(wvice>max(iStruct.bandPass(:)),1,'first');
waves = wvice(k1:k2);

%output values
R = zeros(height(iStruct.Tbl),size(iStruct.bandPass,1));

% cycle through all the bands
T = iStruct.Tbl;
for b=1:size(iStruct.bandPass,1)
    k1 = find(waves<min(iStruct.bandPass(b,:)),1,'last');
    k2 = find(waves>max(iStruct.bandPass(b,:)),1,'first');
    thisWave = waves(k1:k2);
    sol = SolarScale(thisWave,'units',iStruct.waveUnit);
    FS = fit(thisWave,sol,'pchipinterp');
    % cycle through the radii
    integralR = zeros(height(iStruct.Tbl),1);
    for k=1:height(T)
        specArg = buildArgs(iStruct,k,b);
        if isempty(specArg)
            [thisR,rWave] = SnowCloudSpectralReflectance(T.cosZ(k),iStruct.substance,...
                T.radius(k),iStruct.sizeUnit,thisWave,iStruct.waveUnit);
        else
            [thisR,rWave] = SnowCloudSpectralReflectance(T.cosZ(k),iStruct.substance,...
                T.radius(k),iStruct.sizeUnit,thisWave,iStruct.waveUnit,...
                specArg{:});
        end
        F = fit(rWave,thisR,'pchipinterp');
        if iStruct.ignoreSolar
            integralR(k) = integrate(F,iStruct.bandPass(b,2),iStruct.bandPass(b,1))/...
                (iStruct.bandPass(b,2)-iStruct.bandPass(b,1));
        else
            y = F(thisWave).*FS(thisWave);
            Fy = fit(thisWave,y,'pchipinterp');
            integralR(k) = integrate(Fy,iStruct.bandPass(b,2),iStruct.bandPass(b,1))/...
                integrate(FS,iStruct.bandPass(b,2),iStruct.bandPass(b,1));
        end
        
    end
    R(:,b) = integralR;
end
end

function iStruct = parseInput(cosZ,substance,radius,radiusUnit,varargin)
%parse input values

% some parameters
S = SnowCloudLimits;
defaultWE = Inf;
defaultWEunits = 'mm';
defaultR0 = 0;
waveUnit = 'nm';
sizeUnit = 'mum';

p = inputParser;
rangeValidation = @(x) isnumeric(x) && all(x(:)>=0 & x(:)<=1);
positiveValidation = @(x) isnumeric(x) && all(x(:)>0);
bandValidation = @(x) (isrow(x) || iscolumn(x)) &&...
    ((isnumeric(x) && all(x(:)>0)) || iscell(x) || iscategorical(x));
bpValidation = @(x) isnumeric(x) && all(x(:)>=0) &&...
    (isvector(x) || size(x,2)==2);
addRequired(p,'cosZ',rangeValidation)
addRequired(p,'substance',@ischar);
addRequired(p,'radius',@isnumeric)
addRequired(p,'radiusUnit',@ischar)
addParameter(p,'we',defaultWE,positiveValidation)
addParameter(p,'weunits',defaultWEunits,@ischar)
addParameter(p,'waveunit',waveUnit,@ischar)
addParameter(p,'r0',defaultR0,rangeValidation)
addParameter(p,'dust',[],rangeValidation)
addParameter(p,'dustradius',S.defaultDustRadius,positiveValidation)
addParameter(p,'soot',[],rangeValidation)
addParameter(p,'sootradius',S.defaultSootRadius,positiveValidation)
addParameter(p,'cloudradius',S.defaultWaterCloudRadius,positiveValidation)
addParameter(p,'lookup',true,@islogical)
addParameter(p,'waterconc',0,rangeValidation)
addParameter(p,'bandpass',[],bpValidation)
addParameter(p,'sensor','',@ischar)
addParameter(p,'band',[],bandValidation)
addParameter(p,'ignoresolar',false,@islogical)
parse(p,cosZ,substance,radius,radiusUnit,varargin{:})

% check if spectrometer, which will trigger a direct call to
% SnowCloudSpectralReflectance
iStruct.spectrometer = false;
if strcmpi(p.Results.sensor,'aviris-ng')
    iStruct.spectrometer = true;
    T = SensorTable('aviris-ng',p.Results.waveunit);
    iStruct.wavelength = T.CentralWavelength;
    iStruct.waveUnit = p.Results.waveunit;
elseif ~isempty(p.Results.bandpass)
    if isvector(p.Results.bandpass)
        iStruct.spectrometer = true;
        iStruct.wavelength = p.Results.bandpass;
        iStruct.waveUnit = p.Results.waveunit;
    end
end
if iStruct.spectrometer
    return
end

% some logical variables
iStruct.wetSnow = contains(p.Results.substance,'wetsnow','IgnoreCase',true);
iStruct.mixedPhase = (contains(p.Results.substance,'mixedphase','IgnoreCase',true));
iStruct.dustySnow = ~isempty(p.Results.dust);
iStruct.sootySnow = ~isempty(p.Results.soot);
iStruct.deepSnow = all(isinf(p.Results.we(:)));
iStruct.cleanSnowCloud = ~(iStruct.dustySnow || iStruct.sootySnow);

% units
iStruct.sizeUnit = sizeUnit;
iStruct.waveUnit = waveUnit;

% begin transfer to data structure
% sensor/band characteristics
iStruct.bandPass = getBands(p,iStruct);
iStruct.ignoreSolar = p.Results.ignoresolar;

% snow/cloud characteristics
substance = lower(p.Results.substance);
assert(strcmp(substance,'ice') || strcmp(substance,'snow') ||...
    strcmp(substance,'water') || strcmp(substance,'wetsnow') ||...
    strncmp(substance,'mixedphase',10),...
    'substance ''%s'' not recognized',p.Results.substance)
iStruct.substance = substance;

% check for some incompatibilities
if iStruct.wetSnow
    assert(~isempty(p.Results.waterconc) && all(p.Results.waterconc(:)<=0.15),...
        'if ''substance'' is ''wetSnow'', ''waterConc'' must be <= 0.15')
elseif iStruct.mixedPhase
    assert(~isempty(p.Results.waterconc),...
        'if ''substance'' is ''mixedPhaseCloud'', ''waterConc'' must be specified and ''cloudRadius'' should be specified')
end

% check and adjust sizes
[mu0,radius,WE] = checkSizes(p.Results.cosZ,p.Results.radius,...
    p.Results.we);
[iStruct.R0,~] = checkSizes(p.Results.r0,iStruct.bandPass(:,1));

% convert in input are matrices, save to reconvert later
iStruct.origSize = size(mu0);
T = table(mu0(:),convertUnits(radius(:),p.Results.radiusUnit,iStruct.sizeUnit),...
    WE(:),'VariableNames',{'cosZ','radius','WE'});

if ~iStruct.cleanSnowCloud
    if iStruct.dustySnow
        if p.Results.dustradius == S.defaultDustRadius
            dustRadius = convertUnits(p.Results.dustradius,S.unitsSize,...
                iStruct.sizeUnit);
        else
            dustRadius = convertUnits(p.Results.dustradius,...
                p.Results.radiusUnit,iStruct.sizeUnit);
        end
        [dustConc,dustRadius,~,~] =...
            checkSizes(p.Results.dust,dustRadius,radius,WE);
        newT = table(dustConc,dustRadius);
        T = [T newT];
    end
    if iStruct.sootySnow
        sootConc = p.Results.soot;
        if p.Results.sootradius == S.defaultSootRadius
            sootRadius = convertUnits(p.Results.sootradius,S.unitsSize,'mum');
        else
            sootRadius = convertUnits(p.Results.sootradius,...
                p.Results.radiusUnit,'mum');
        end
        [sootConc,sootRadius,~,~] =...
            checkSizes(sootConc,sootRadius,radius,WE);
        newT = table(sootConc,sootRadius);
        T = [T newT];
    end
end

% mixed phase or wet snow
if iStruct.wetSnow || iStruct.mixedPhase
    waterConc = p.Results.waterconc;
    if iStruct.mixedPhase
        if p.Results.cloudradius==S.defaultWaterCloudRadius
            cloudRadius = convertUnits(S.defaultWaterCloudRadius,S.unitsSize,...
                iStruct.sizeUnit);
        else
            cloudRadius = convertUnits(p.Results.cloudradius,...
                p.Results.radiusUnit,iStruct.sizeUnit);
        end
        [waterConc,cloudRadius,~] = checkSizes(waterConc,cloudRadius,radius);
        newT = table(waterConc,cloudRadius);
    else
        [waterConc,~] = checkSizes(waterConc,radius);
        newT = table(waterConc);
    end
    T = [T newT];
end

% eliminate duplicates
T = unique(T);
iStruct.Tbl = T;

% use Mie Lookup tables
iStruct.lookup = p.Results.lookup;
end

function bandPass = getBands(p,iStruct)
% check consistence of band designations
assert(xor(isempty(p.Results.bandpass),isempty(p.Results.sensor)),...
    'either ''bandPass'' or ''sensor''/''bands'' must be specified')
if ~isempty(p.Results.bandpass)
    bandPass = convertUnits(p.Results.bandpass,p.Results.waveunit,iStruct.waveUnit);
else
    T = SensorTable(p.Results.sensor,iStruct.waveUnit);
    if isempty(p.Results.band)
        % all bands
        bandPass = [T.LowerWavelength T.UpperWavelength];
    else
        x = p.Results.band;
        if isnumeric(x)
            band = categorical(x);
        elseif iscategorical(x)
            band = x;
        else % cell
            band = categorical(x);
        end
        bandPass = zeros(length(x),2);
        for k=1:length(band)
            b = find(T.Band==band(k));
            if isempty(b)
                warning(['band ' band(k) ' not found'])
            end
            bandPass(k,:) = [T.LowerWavelength(b) T.UpperWavelength(b)];
        end
    end
end
% make sure bigger number is in column 2
bandPass = sort(bandPass,2);
end

function argc = buildArgs(iStruct,k,b)
% build argument list for passing to SnowCloudSpectralReflectance
% k is index into mu0, b is index into bands

argc = {''};
T = iStruct.Tbl;
if iStruct.wetSnow || iStruct.mixedPhase
    argc = cat(2,argc,'waterConc',T.waterConc(k));
    if iStruct.mixedPhase && isfield(T,'cloudRadius')
        argc = cat(2,argc,'cloudRadius',T.cloudRadius(k));
    end
end
if ~iStruct.cleanSnowCloud
    if iStruct.dustySnow
        argc = cat(2,argc,'dust',T.dustConc(k));
        if isfield(T,'dustRadius')
            argc = cat(2,argc,'dustRadius',T.dustRadius(k));
        end
    end
    if iStruct.sootySnow
        argc = cat(2,argc,'soot',T.sootConc(k));
        if isfield(T,'sootRadius')
            argc = cat(2,argc,'sootRadius',T.sootRadius(k));
        end
    end
end
if ~iStruct.deepSnow
    argc = cat(2,argc,'WE',T.WE(k));
    if isfield(iStruct,'R0')
        argc = cat(2,argc,'R0',iStruct.R0(b));
    end
end
% eliminate the leading empty cell
argc = argc(2:end);
end