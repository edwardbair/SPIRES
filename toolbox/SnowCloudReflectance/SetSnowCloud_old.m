function fscript = SetSnowCloud(callingFunction,varargin)
% usage 1: prescription=SetSnowCloud(mfilename,substance,prop1,val1,pro2,val2,...)
% usage 2: prescription=SetSnowCloud(mfilename,substance,fscript,prop1,val1,pro2,val2,...)
%SetSnowCloud: defines the prescription for a reflectance model based on
%snow or cloud properties
% A set of property/value pairs are parsed into a structure to use as the
% prescription for all the functions that calculate reflectance,
% emissivity, and transmittance
%
% (inspired by John D'Errico's slmset for the SLM toolbox, available on
% the MATLAB File Exchange)
%
%First argument, required, is the calling function, usually by specifying
%mfilename
%
%substance and fscript argument must come before the name-value pairs if
%either or both are entered, in either order
%   substance, either 'snow', 'iceCloud', 'waterCloud', or 'mixedCloud'
%       (but any unambiguous abbreviation beginning with first letter works)
%   fscript, if you want to modify an existing prescription
%
%%remaining arguments
%   prop/val as a set of name-value pairs (detailed descriptions below)
%       property names are case-insensitive and can be shortened as long as
%       the short name is unambiguous
%   'cosZ' - cosine of illumination angle on flat surface, scalar, vector, or matrix
%   'muS' - cosine of illumination angle on slope
%   'cosS' - cosine of slope
%   'radius', effective optical radius of snow grain, cloud ice crystal, or
%       cloud water droplet
%   alternatively, enter 'ssa' and specific surface area, in m^2/kg, if
%       substance is 'snow', 'iceCloud', or 'mixedCloud'
%%
% Wavelength properties, either 'wavelength', 'bandPass', or 'sensor'/'band'
%   must be specified, but not both
% 'wavelength' – vector if sensor is a spectrometer
% 'bandPass' Nx2 matrix if sensor is multispectral (use [0.28 4] or [280 4000]
%   to get full-spectrum albedo)
% 'waveUnit' – units for wavelength, no default (to prevent errors)
% 'sensor' – instead of 'wavelength' or 'bandPass', can specify spectrometer
%   or multispectral sensor, anything in the SensorTable.m function
%   (if calling function is 'spectral', returns CentralWavelength, a
%   vector, or if calling function is 'intg', returns Nx2 matrix)
% 'bands' – bands of the sensor, either numeric vector, cell vector, or
%   categorical vector of bands,or, if omitted, all bands for that sensor
% 'ignoreSolar' – if false (default), solar radiation accounted for unless
%   outside range of SolarScale.m
%   if true, ignores solar radiation and just provides band-average
%   reflectivity (this is needed to calculate emissivity)
%   set to false automatically if sensor is a spectrometer
%%
% Properties applicable to either snow or cloud
% 'sizeUnit' – units for optically equivalent radius of snow grains (default
%   'mum', but any common metric length unit will be converted)
% 'dust' – mass fraction
% 'dustRadius' – same units as for optically equivalent snow grain radius
% 'soot' – mass fraction
% 'sootRadius' – same units as for optically equivalent snow grain radius
% 'waterEquivalent' – water equivalent, Inf if not specified but must be specified for cloud
% 'weUnit' – unit for measuring WE, default 'mm'
% 'R0' – reflectance of surface under cloud or snow, or, if 'WE' is not
%   specified (or specified as Inf), then treat as a fractional snow mixed
%   with dirt and/or vegetation.
%   Specify as a scalar or, more commonly, as a table with columns wavelength
%   and reflectance, with units for wavelength specified in
%   Tbl.Properties.VariableUnits (assume same as 'waveUnit' if not specified)
%   For fractional snow, reflectance can have one or two columns depending
%   on whether pixel has one or two non-snow endmembers
% 'temperature' - degrees K, physical temperature of snow or cloud, used
%   only in calculations involving wavelength-integrated emissivity
%%
% Properties applicable only to snow
% 'wetness' – water mass fraction (0 to 0.2)
% 'fractionalCoverage' – if 'waterEquivalent' is not specified (or is Inf)
%   and R0 is a scalar or a vector, then a 2-element vector [fSCA fOther]
%   that sums to 1.0, or if
%   R0 is a matrix, then a 3-element vector [fSCA fMem1 fMem2] that also
%   sums to 1.0 (Mem1 and Mem2 might be soil and vegetation, for example)
% alternatively, 'fSCA' can be specified instead, in which case
%   fractionalCoverage is set to [fSCA 1-fSCA]
%%
% Properties applicable only to mixed clouds
% In this case, the required radius argument is the size of the ice crystals
% 'wetness' – if 'substance' is 'mixed', water mass fraction (0 to 1)
% 'waterRadius' must be specified if 'mixed' (the required radius argument is
% the size of the ice crystals)
%%
% Properties about the radiative transfer calculations
% 'lookup' – Use lookup tables to calculate Mie variables, default true
% 'method' - (default is delta-Eddington, can specify by abbreviations)
%   'meador' or 'hybrid' - for Meader-Weaver hybrid
%   'delta' or 'eddington' - for delta-Eddington (like Wiscombe & Warren)
%   'disort' (discrete ordinates)'
% (generally any 3- or more-letter abbreviation works)
%% Properties about the inversion
% 'solutionmethod' - default is 'lsqnonlin', alternative is 'spectralangle'
%   (default of 3 or more letters works)

%%
persistent lastSensor lastWavelength lastBandPass
narginchk(1,Inf)
nargoutchk(0,1)
S = SnowCloudLimits;
strings = {'snow','ice','iceCloud','water','waterCloud','mixed',...
    'mixedCloud','mixedPhase','mixedPhaseCloud'};

% process inputs (complicated, lots of options)
p = inputParser;

% get calling function to later see if it's an inverse problem
addRequired(p,'callingFunction',@ischar)

%% are substance and/or an initial prescription struct provided?
chooseDefaultP = true;
chooseSubstance = false;
rmk = false(1,min(2,length(varargin)));
for k=1:min(2,length(varargin))
    if isstruct(varargin{k})
        fscript = varargin{k};
        chooseDefaultP = false;
        rmk(k) = true;
    elseif ischar(varargin{k}) && any(strcmpi(varargin{k},strings))
        substance = varargin{k};
        rmk(k) = true;
        chooseSubstance = true;
    end
end

% get rid of the variable arguments if not defaults
getRid = sum(rmk);
if getRid>0
    varargin(1:getRid) = [];
end
if isempty(varargin)
    return
end

%% all substances
snow = categorical({'snow'});
iceCloud = categorical({'iceCloud'});
waterCloud = categorical({'waterCloud'});
mixedCloud = categorical({'mixedCloud'});
notSet = categorical({'notSet'}); %#ok<NASGU>

%% set the substance, which might have changed from previously existing prescription
if chooseSubstance
    matchstr = validatestring(substance,strings);
    switch matchstr
        case 'snow'
            fscript.substance = snow;
            fscript.chooseWater = false;
        case {'ice','iceCloud'}
            fscript.substance = iceCloud;
            fscript.chooseWater = false;
        case {'water','waterCloud'}
            fscript.substance = waterCloud;
            fscript.chooseWater = true;
        case {'mixed','mixedCloud','mixedPhase','mixedPhaseCloud'}
            fscript.substance = mixedCloud;
            fscript.chooseWater = false;
        otherwise
            error('No match for ''%s''',matchstr)
    end
end

% select the default prescription if not provided of if substance has changed
if chooseDefaultP
    % set defaults -- wavelength or sensor
    fscript.wavelength = [];
    fscript.bandPass = [];
    fscript.waveUnit = '';
    fscript.sensor = [];
    fscript.bands = [];
    fscript.ignoreSolar = false;
    fscript.cosZ = [];
    fscript.muS = [];
    fscript.cosS = [];
    % used only in inverse problem, from eqn Reff = R*cosS*muS/mu0;
    % corrFactor = cosS*muS
    fscript.corrFactor = [];
    
    % defaults - snow or cloud
    fscript.sizeUnit = 'mum';
    if fscript.substance==snow
        fscript.WE = Inf;
    else
        fscript.WE = [];
    end
    fscript.weUnit = 'mm';
    fscript.R0 = 0;
    fscript.waterRadius = S.defaultWaterCloudRadius;
    if fscript.substance==snow
        fscript.iceRadius = S.defaultSnowRadius;
    else
        fscript.iceRadius = S.defaultIceCloudRadius;
    end
    fscript.SSA = radius2SSA(fscript.iceRadius,fscript.sizeUnit);
    fscript.wetness = [];
    fscript.fractionalCoverage = [];
    fscript.temperature = 273.16;
    
    % defaults -- dust and soot
    fscript.dust = [];
    fscript.dustRadius = S.defaultDustRadius;
    fscript.soot = [];
    fscript.sootRadius = S.defaultSootRadius;
    
    % defaults -- radiative transfer
    fscript.lookup = true;
    fscript.method = 'delta-Eddington';
    
    % defaults -- solution method for inversion
    fscript.doInverse = false;
    fscript.solutionMethod = 'lsqnonlin';
end

%% parse input and check sizes
bandValidation = @(x) ((isnumeric(x) && all(x(:)>0)) || iscell(x) || iscategorical(x));
bpValidation = @(x) isnumeric(x) && all(x(:)>=0);
r0Validation = @(x) (isnumeric(x) && isscalar(x) && x>=0 && x<=1) ||...
    (istable(x) &&...
    ~isempty(contains(x.Properties.VariableNames,'wavelength','IgnoreCase',true)) &&...
    ~isempty(contains(x.Properties.VariableNames,'reflectance','IgnoreCase',true))) ||...
    (isnumeric(x) && all(x(:)>=0) && all(x(:)<=1));
nonNegValidation = @(x) isnumeric(x) && all(x(:)>=0) && all(x(:)<=1);
positiveValidation = @(x) isnumeric(x) && all(x(:)>0);

%illumination angle(s)
addParameter(p,'cosZ',fscript.cosZ,nonNegValidation)
addParameter(p,'muS',fscript.muS,nonNegValidation)
addParameter(p,'cosS',fscript.cosS,nonNegValidation)
addParameter(p,'corrfactor',fscript.corrFactor,@isnumeric)

% wavelength or sensor (can specify just one, not both)
% radius or SSA (can specify just one, not both)
iswave = false;
issens = false;
isbandp = false;
isradius = false;
isSSA = false;
fscript.spectrometer = contains(callingFunction,'spectral','IgnoreCase',true);
for k=1:length(varargin)
    if ischar(varargin{k})
        if ~iswave
            iswave = contains(varargin{k},'wavel','IgnoreCase',true);
            if iswave % reading new wavelength, so set default to empty
                fscript.wavelength = [];
            end
        end
        if ~issens
            issens = contains(varargin{k},'sens','IgnoreCase',true);
            if issens % reading new sensor, so set defaults to empty
                fscript.sensor = [];
                fscript.bands = [];
            end
        end
        if ~isbandp
            isbandp = contains(varargin{k},'bandp','IgnoreCase',true);
            if isbandp % reading new bandPass, so set defaults to empty
                fscript.bandPass = [];
            end
        end
        if ~isradius
            isradius = strcmpi(varargin{k},'radius');
            if isradius % reading new radius, so set defaults to empty
                switch fscript.substance
                    case {snow,iceCloud,mixedCloud}
                        fscript.iceRadius = [];
                        fscript.SSA = [];
                    otherwise
                        fscript.waterCloudRadius = [];
                end
            end
        end
        if ~isSSA
            isSSA = strcmpi(varargin{k},'ssa');
            if isSSA % reading new SSA, so set defaults to empty
                fscript.SSA = [];
                fscript.iceRadius = [];
            end
        end
    end
end
assert(sum([iswave issens isbandp])<=1,...
    '''sensor'' or ''wavelength'' or ''bandPass'' can be specified, but not more than one')
addParameter(p,validatestring('wavelength',{'wavel','wavelength'}),...
    fscript.wavelength,bpValidation)
addParameter(p,validatestring('bandpass',{'bandp','bandpass'}),...
    fscript.bandPass,bpValidation)
addParameter(p,validatestring('waveunit',{'waveu','waveunit','lambda'}),...
    fscript.waveUnit,@ischar)
addParameter(p,validatestring('sensor',{'sens','sensor'}),...
    fscript.sensor,@(x) ischar(x) || iscategorical(x))
addParameter(p,validatestring('bands',{'band','bands'}),...
    fscript.bands,bandValidation)
addParameter(p,validatestring('ignoresolar',{'ign','ignore','ignoresolar'}),...
    fscript.ignoreSolar,@islogical)

% snow or cloud
addParameter(p,validatestring('sizeunit',{'size','sizeunit'}),...
    fscript.sizeUnit,@ischar)
addParameter(p,'dust',fscript.dust,@(x) isscalar(x) && x>=0 && x<1);
addParameter(p,validatestring('dustradius',{'dustr','dustradius','dradius'}),...
    fscript.dustRadius,@(x) isscalar(x) && x>0)
addParameter(p,'soot',fscript.soot,@(x) isscalar(x) && x>=0 && x<1);
addParameter(p,validatestring('sootradius',{'sootr','sootradius','sradius'}),...
    fscript.sootRadius,@(x) isscalar(x) && x>0)
addParameter(p,validatestring('waterequivalent',{'waterequiv','waterequivalent','wequiv'}),...
    fscript.WE,@(x) isscalar(x) && x>0)
addParameter(p,validatestring('weunit',{'weu','weunit'}),...
    fscript.weUnit,@ischar)
addParameter(p,'r0',fscript.R0,r0Validation)
addParameter(p,validatestring('temperature',{'tem','temp','temperature'}),...
    fscript.temperature,@(x) isscalar(x) && x>0)

% some parameters to load, even if not used
addParameter(p,validatestring('fractionalcoverage',{'frac','fraccover',...
    'fractionalcoverage','fraction'}),fscript.fractionalCoverage,...
    @(x) isnumeric(x) && (length(x)==2 || length(x)==3))
addParameter(p,'fsca',[],@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
addParameter(p,validatestring('waterradius',{'water','waterr','waterradius'}),...
    S.defaultWaterCloudRadius,@(x) isscalar(x) && x>0);

% defaults - snow
assert(~(isradius && isSSA),'''radius'' or ''SSA'' can be specified, but not both')
if fscript.substance==snow
    addParameter(p,'wetness',fscript.wetness,@(x) isscalar(x) &&...
        x>=S.wetSnow(1) && x<=S.wetSnow(2))
    addParameter(p,'radius',fscript.iceRadius,positiveValidation)
    addParameter(p,'ssa',fscript.SSA,positiveValidation)
    % defaults - cloud
elseif fscript.substance==iceCloud
    addParameter(p,'radius',fscript.iceRadius,positiveValidation)
    addParameter(p,'ssa',fscript.SSA,positiveValidation)
elseif fscript.substance==waterCloud
    addParameter(p,'radius',fscript.waterRadius,positiveValidation)
elseif fscript.substance==mixedCloud
    addParameter(p,'wetness',fscript.wetness,@(x) isscalar(x) && x>=0 && x<=1)
    addParameter(p,'radius',fscript.iceRadius,positiveValidation)
end

% radiative-transfer
addParameter(p,validatestring('lookup',{'look','lookup'}),...
    fscript.lookup,@islogical)
addParameter(p,validatestring('method',{'meth','method'}),...
    fscript.method,@ischar)

% inversion
addParameter(p,validatestring('solutionmethod',{'sol','soln','solu','solnmethod','solutionmethod'}),...
    fscript.solutionMethod,@ischar)

% insert into structure
parse(p,callingFunction,varargin{:})

% inverse problem?
callingFunction = p.Results.callingFunction;
fscript.doInverse = contains(callingFunction,'inver','IgnoreCase',true);

% check for sizeUnit, will be converted to 'mum'
if ~(strcmpi(p.Results.sizeunit,'um') || strcmpi(p.Results.sizeunit,'mum'))
    warning('sizes of scatterers will be converted to micrometers')
end


% warning message about unnecessary (thereby ignored) entries
if contains(char(fscript.substance),'cloud','IgnoreCase',true) &&...
        (~isempty(p.Results.fractionalcoverage) || ~isempty(p.Results.fsca))
    warning('''fractionalCoverage'' or ''fSCA'' not appropriate for clouds, ignored')
end

% check whether to use sensor, wavelengths, or bandPass
newWavelength = ~issens && ~isbandp &&...
    (iswave || isempty(lastWavelength) || ~isequal(p.Results.wavelength,lastWavelength));
newSensor = ~iswave && ~isbandp &&...
    (issens || isempty(lastSensor) || fscript.sensor~=lastSensor);
newBandPass = ~iswave && ~issens &&...
    (isbandp || isempty(lastBandPass) || ~isequal(p.Results.bandpass,lastBandPass));
if newWavelength
    fscript.sensor = [];
    fscript.bands = [];
    fscript.bandPass = [];
    fscript.wavelength = p.Results.wavelength;
elseif newBandPass
    fscript.sensor = [];
    fscript.bands = [];
    fscript.bandPass = p.Results.bandpass;
    lastBandPass = fscript.bandPass;
elseif newSensor
    if iscategorical(p.Results.sensor)
        fscript.sensor = p.Results.sensor;
    else
        fscript.sensor = categorical({p.Results.sensor});
    end
    if iscategorical(p.Results.bands)
        fscript.bands = p.Results.bands;
    else
        fscript.bands = categorical(p.Results.bands);
    end
end

fscript.ignoreSolar = p.Results.ignoresolar;
fscript.waveUnit = p.Results.waveunit;

% if sensor specified, convert to wavelengths or bandPass
if ~isempty(fscript.sensor)
    Tbl = SensorTable(fscript.sensor,fscript.waveUnit);
    fscript.sensor = unique(Tbl.Sensor);
    % spectrometer
    if contains(callingFunction,'spectral','IgnoreCase',true)
        fscript.wavelength = Tbl.CentralWavelength;
        fscript.bands = Tbl.Band';
        fscript.spectrometer = true;
        fscript.multispectral = false;
        fscript.bandPass = [Tbl.LowerWavelength Tbl.UpperWavelength];
    elseif contains(callingFunction,'intg','IgnoreCase',true) % multispectral?
        % make sure we're not using the bands from the previous sensor
        if isempty(p.Results.bands)
            fscript.bandPass = [Tbl.LowerWavelength Tbl.UpperWavelength];
            fscript.bands = Tbl.Band';
        else
            x = p.Results.bands;
            if isnumeric(x)
                fscript.bands = categorical(x);
            elseif iscategorical(x)
                fscript.bands = x;
            else % cell
                fscript.bands = categorical(x);
            end
            bandPass = zeros(length(x),2);
            wv = zeros(length(x),1);
            for k=1:length(fscript.bands)
                b = find(Tbl.Band==fscript.bands(k));
                if isempty(b)
                    warning(['band ' fscript.bands(k) ' not found'])
                end
                bandPass(k,:) = [Tbl.LowerWavelength(b) Tbl.UpperWavelength(b)];
                wv(k) = Tbl.CentralWavelength(b);
            end
            fscript.bandPass = bandPass;
            fscript.wavelength = wv;
        end
        fscript.spectrometer = false;
        fscript.multispectral = true;
    end
elseif ~isempty(p.Results.bandpass)
    fscript.bandPass = p.Results.bandpass;
    assert(size(fscript.bandPass,2)==2,...
        '''bandPass'' must be a matrix with 2 columns')
    fscript.spectrometer = false;
    fscript.multispectral = true;
elseif ~isempty(p.Results.wavelength)
    fscript.wavelength = p.Results.wavelength;
    assert(isvector(fscript.wavelength),...
        '''wavelength'' must be a vector');
    fscript.spectrometer = true;
    fscript.multispectral = false;
end

% if spectrum, radius and cosZ and wavelength must be same size, but if bandpass
% sizes of cosZ and radius can be arbitrary but are converted to 1D unique
% vectors
if isfield(fscript,'spectrometer') && fscript.spectrometer
    % conversion actually not necessary, but in place for future mods
    if isempty(p.Results.cosZ)
        if isSSA
            [radius,fscript.wavelength] =...
                checkSizes(SSA2radius(p.Results.ssa,S.unitsSize),...
                fscript.wavelength);
        else
            [radius,fscript.wavelength] =...
                checkSizes(convertLengthUnits(p.Results.radius,...
                p.Results.sizeunit,S.unitsSize),fscript.wavelength);
        end
    else
        if isSSA
            [fscript.cosZ,radius,fscript.wavelength] = checkSizes(p.Results.cosZ,...
                SSA2radius(p.Results.ssa,S.unitsSize),...
                fscript.wavelength);
        else
            [fscript.cosZ,radius,fscript.wavelength] =...
                checkSizes(p.Results.cosZ,...
                convertLengthUnits(p.Results.radius,p.Results.sizeunit,S.unitsSize),...
                fscript.wavelength);
            
        end
        fscript.cosZ = unique(fscript.cosZ);
    end
else
    fscript.cosZ = unique(p.Results.cosZ(:));
    % conversion actually not necessary, but in place for future mods
    if isSSA
        radius = SSA2radius(unique(p.Results.ssa(:)),S.unitsSize);
    else
        radius = convertLengthUnits(unique(p.Results.radius(:)),...
            p.Results.sizeunit,S.unitsSize);
    end
end
if isempty(p.Results.muS) && isempty(p.Results.cosS)
    fscript.ignoreSlope = true;
else
    assert(fscript.substance==snow,...
        'inappropriate to specify ''muS'' and/or ''cosS'' unless substance is snow')
    assert(~isempty(fscript.cosZ),...
        'if ''muS'' and/or ''cosS'' is specified, ''cosZ'' must also')
    [fscript.muS,fscript.cosS,~] = checkSizes(p.Results.muS,...
        p.Results.cosS,fscript.cosZ);
    fscript.ignoreSlope = false;
end
switch fscript.substance
    case snow
        fscript.iceRadius = radius;
    case iceCloud
        fscript.iceRadius = radius;
    case waterCloud
        fscript.waterRadius = radius;
    case mixedCloud
        fscript.iceRadius = radius;
        fscript.waterRadius =...
            convertLengthUnits(p.Results.waterradius,p.Results.sizeunit,S.unitsSize);
    otherwise
        error(['fscript.substance=' fscript.substance ' not recognized']) % shouldn't reach
end

% snow or cloud
R0 = p.Results.r0;
if istable(R0)
    assert(all(R0.wavelength>0) && all(R0.reflectance>=0) &&...
        all(R0.reflectance<=1),'R0 table bad wavelength or reflectance')
elseif isnumeric(R0) && ~isscalar(R0)
    assert(isequal(size(R0,1),size(fscript.wavelength,1)) ||...
        isequal(length(R0),length(fscript.wavelength)),...
        'R0 if vector must be same size as wavelength vector or matrix (dimension 1)')
end
fscript.R0 = R0;
fscript.sizeUnit = S.unitsSize;
if isempty(p.Results.dust)
    fscript.dustySnowCloud = false;
else
    fscript.dustySnowCloud = true;
    fscript.dust = p.Results.dust;
    if ~strcmpi(S.unitsSize,fscript.sizeUnit) &&...
            p.Results.dustradius==S.defaultDustRadius
        % conversion actually not necessary, but in place for future mods
        fscript.dustRadius = convertLengthUnits(p.Results.dustradius,S.unitsSize,fscript.sizeUnit);
    else
        fscript.dustRadius = p.Results.dustradius;
    end
end
if isempty(p.Results.soot)
    fscript.sootySnowCloud = false;
else
    fscript.sootySnowCloud = true;
    fscript.soot = p.Results.soot;
    if ~strcmpi(S.unitsSize,fscript.sizeUnit) &&...
            p.Results.sootradius==S.defaultSootRadius
        % conversion actually not necessary, but in place for future mods
        fscript.sootRadius = convertLengthUnits(p.Results.sootradius,S.unitsSize,fscript.sizeUnit);
    else
        fscript.sootRadius = p.Results.sootradius;
    end
end

% snow
if fscript.substance==snow
    if isempty(p.Results.wetness)
        fscript.wetSnow = false;
    else
        fscript.wetSnow = true;
        if ~isempty(p.Results.wetness)
            fscript.wetness = p.Results.wetness;
        end
    end
    if isinf(p.Results.waterequivalent)
        fscript.WE = p.Results.waterequivalent;
        fscript.deepSnow = true;
    else
        fscript.deepSnow = false;
        fscript.WE = convertLengthUnits(p.Results.waterequivalent,p.Results.weunit,'mm');
        fscript.weUnit = 'mm';
    end
    if ~isempty(p.Results.fsca)
        fscript.fractionalCoverage = [p.Results.fsca 1-p.Results.fsca];
        if ~isempty(p.Results.fractionalcoverage)
            warning('fractionalCoverage overwritten by fSCA')
        end
    else
        fscript.fractionalCoverage = p.Results.fractionalcoverage/...
            sum(p.Results.fractionalcoverage);
    end
elseif contains(char(fscript.substance),'cloud','IgnoreCase',true)
    if ~fscript.doInverse
        assert(~isempty(p.Results.waterequivalent) &&...
            ~isinf(p.Results.waterequivalent),...
            'if substance is some sort of cloud, you must specify a non-infinite ''waterEquivalent''')
    end
    fscript.WE = convertLengthUnits(p.Results.waterequivalent,p.Results.weunit,'mm');
    fscript.weUnit = 'mm';
    if fscript.substance==mixedCloud && ~fscript.doInverse
        assert(~isempty(p.Results.wetness),...
            'for ''mixedCloud'' you must also specify ''wetness''')
        fscript.wetness = p.Results.wetness;
    end
else
    error(['substance=' fscript.substance ' not recognized']); % shouldn't reach
end

% radiative transfer
method = p.Results.method;
if strncmpi(method,'two',3) || strncmpi(method,'mea',3) || strncmpi(method,'hyb',3)
    fscript.method = 'Meador-Weaver hybrid twostream';
elseif strncmpi(method,'del',3) || strncmpi(method,'edd',3)
    fscript.method = 'delta-Eddington twostream';
elseif strncmpi(method,'dis',3)
    fscript.method = 'disort';
else
    error('''method'' = ''%s'' unrecognized',p.Results.method)
end

% solution method if inversion
if strncmpi(p.Results.solutionmethod,'lsq',3)
    fscript.solutionMethod = 'lsqnonlin';
elseif strncmpi(p.Results.solutionmethod,'spe',3)
    fscript.solutionMethod = 'spectralAngle';
else
    error('''solutionMethod'' ''%s'' not recognized',p.Results.solutionmethod)
end

%Mie lookup
fscript.lookup = p.Results.lookup;

% add SSA or radius to the output
if isfield(fscript,'iceRadius')
    if ~isempty(fscript.iceRadius)
        fscript.SSA = radius2SSA(fscript.iceRadius,fscript.sizeUnit);
    end
elseif isfield(fscript,'SSA')
    if ~isempty(fscript.SSA)
        fscript.iceRadius = SSA2radius(fscript.SSA,fscript.sizeUnit);
    end
    
end

lastSensor = fscript.sensor;
lastWavelength = fscript.wavelength;
end