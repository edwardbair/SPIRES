function fscript = SetSpectral(callingFunction,varargin)
% usage 1: prescription=SetSpectral(mfilename,substance,prop1,val1,pro2,val2,...)
% usage 2: prescription=SetSpectral(mfilename,substance,fscript,prop1,val1,pro2,val2,...)
%SetSpectral: defines the prescription for a reflectance model based on
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
%%  illumination geometry
%   'cosZ' - cosine of illumination angle on surface, scalar, vector, or matrix
%       (set to empty to calculate diffuse reflectance)
%   'muS' - cosine of illumination angle on slope (default same as cosZ)
%   'viewF' - view factor (default 1)
%% scatterer
%   'radius', effective optical radius of snow grain, cloud ice crystal, or
%       cloud water droplet
%   alternatively, enter 'ssa' and specific surface area, in m^2/kg, if
%       substance is 'snow', 'iceCloud', or 'mixedCloud'
%%
% Wavelength properties, either 'wavelength', 'bandPass', or 'sensor'/'band'
%   must be specified, but just one
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
% 'ignoreSolar' – if false (default), solar radiation accounted for
%   if true, ignores solar radiation and just provides band-average
%   reflectivity (this is needed to calculate emissivity)
%   set to false automatically if sensor is a spectrometer
%%
% Properties applicable to either snow or cloud
% 'sizeUnit' – units for optically equivalent radius of snow grains (default
%   'mum', none others accepted -- might fix sometime)
% 'LAP' - light-absorbing particles, followed either by "dust', 'soot', or
%   'organic', by a cell vector if multiple LAPs are present (i.e. {'dust','soot'},
%   or optionally by a table with 2 columns: wavelength (same units as 'waveUnit')
%   and RefractiveIndex as a complex number with imaginary part positive
% 'Lfraction' – mass fraction of light-absorbing particles, as a vector if
%   multiple LAPS
% 'Lradius' – in same units as for optically equivalent snow grain radius,
%   as a vector if multiple LAPs (can use default if 'dust', 'soot', or 'organic')
% 'waterEquivalent' – water equivalent of snow or cloud in mm (same as kg/m2),
%   Inf if not specified but must be specified for cloud
% 'R0' – reflectance of surface under cloud or snow
%   Specify as a scalar or, more commonly, as a table with columns wavelength
%   and reflectance, with units for wavelength specified in
%   Tbl.Properties.VariableUnits (must specify, will convert if different
%   than waveUnit), or as a griddedInterpolant or cfit object that take
%   wavelength as the argument (in this case code does not check in its
%   wavelength has the correct units)
% 'temperature' - degrees K, physical temperature of snow or cloud, used
%   only in calculations involving wavelength-integrated emissivity
% 'elevation' - in m
%%
% Properties applicable only to snow
% 'wetness' – water mass fraction (0 to 0.2)
%%
% Properties applicable only to snow or mixed clouds
% In this case, the required radius argument is the size of the ice crystals
% 'wetness' – if 'substance' is 'mixed', water mass fraction (0 to 1)
% 'waterRadius' - default 5 um, can be specified if 'mixed' (the required
%   radius argument is the size of the ice crystals)
%%
% Properties about the radiative transfer calculations
% 'icer' - which refractive index data, 'icewb' for Warren & Brandt 2008,
%   'icep' for Picard modification to Warren & Brandt 2008, 'icew' for Warren
%   1984 (default is 'icep')
% 'lookup' – Use lookup tables to calculate Mie variables, default false
% 'method' - (default is delta-Eddington, can specify by abbreviations)
%   'meador' or 'hybrid' - for Meader-Weaver hybrid
%   'delta' or 'eddington' - for delta-Eddington (like Wiscombe & Warren)
%   'disort' (discrete ordinates, not yet implemented)'
% (generally any 3- or more-letter abbreviation works)

%%
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
    fscript.waveUnit = '';
    fscript.sensor = [];
    fscript.cosZ = [];
    fscript.muS = [];
    fscript.viewF = 1;
    
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
    if fscript.substance==waterCloud
        fscript.wetness = [];
    else
        fscript.wetness = 0;
    end
    fscript.temperature = 273.16;
    fscript.icer = 'icep';
    
    % defaults -- light-absorbing particles
    fscript.includeLAP = false;
    fscript.LAP = '';
    fscript.LAPfraction = [];
    fscript.LAPradius = [];
    
    % defaults -- radiative transfer
    fscript.lookupMie = false;
    fscript.method = 'delta-Eddington';
end

%% parse input and check sizes
r0Validation = @(x) (isnumeric(x) && isscalar(x) && x>=0 && x<1) ||...
    (isobject(x) && (isa(x,'griddedInterpolant') || isa(x,'cfit'))) ||...
    (istable(x) &&...
    ~isempty(any(contains(x.Properties.VariableNames,'wavelength','IgnoreCase',true))) &&...
    ~isempty(any(contains(x.Properties.VariableNames,'reflectance','IgnoreCase',true))));
nonNegValidation = @(x) isnumeric(x) && all(x(:)>=0) && all(x(:)<=1);
positiveValidation = @(x) isnumeric(x) && all(x(:)>0);

%illumination angle(s), elevation
addParameter(p,'cosZ',fscript.cosZ,nonNegValidation)
addParameter(p,'muS',fscript.muS,nonNegValidation)
addParameter(p,'viewF',fscript.viewF,nonNegValidation)
addParameter(p,validatestring('elevation',{'elev','elevation'}),3000,...
    @(x) isnumeric(x) && all(x(:)>-500) && all(x(:)<9000))

% wavelength or sensor (can specify just one, not both)
% radius or SSA (can specify just one, not both)
iswave = false;
issens = false;
isradius = false;
isSSA = false;
isR0 = false;

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
        if ~isR0
            isR0 = strcmpi(varargin{k},'R0');
            if isR0
                fscript.R0 = [];
            end
        end
    end
end
assert(sum([iswave issens])<=1,...
    '''sensor'' or ''wavelength'' can be specified, but not more than one')
addParameter(p,validatestring('wavelength',{'wavel','wavelength'}),...
    fscript.wavelength,positiveValidation)
addParameter(p,validatestring('waveunit',{'waveu','waveunit','lambda'}),...
    fscript.waveUnit,@ischar)
addParameter(p,validatestring('sensor',{'sens','sensor'}),...
    fscript.sensor,@(x) ischar(x) || iscategorical(x))
addParameter(p,'bands',[],@(x) ischar(x) || iscategorical(x) || isnumeric(x))
addParameter(p,validatestring('bandpass',{'bandp','bandpass'}),[],...
    @(x) isnumeric(x) && size(x,2)==2)

% snow or cloud
addParameter(p,validatestring('sizeunit',{'size','sizeunit'}),...
    fscript.sizeUnit,@ischar)
addParameter(p,'lap',fscript.LAP,@(x) ischar(x) || iscell(x) || istable(x));
addParameter(p,validatestring('lfraction',{'lfr','lfrac','lfraction','lapfr'}),...
    fscript.LAPfraction,@(x) isnumeric(x) && all(x>=0) && all(x<1))
addParameter(p,validatestring('lradius',{'lr','lrad','lradius','laprad'}),...
    fscript.LAPradius,@(x) isnumeric(x) && all(x>0))
addParameter(p,validatestring('waterequivalent',{'weq','wequiv','waterequiv','waterequivalent'}),...
    fscript.WE,@(x) isscalar(x) && x>0)
addParameter(p,'r0',fscript.R0,r0Validation)
addParameter(p,validatestring('temperature',{'tem','temp','temperature'}),...
    fscript.temperature,@(x) isscalar(x) && x>0)

% some parameters to load, even if not used
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
    fscript.lookupMie,@(x) islogical(x) || isnumeric(x))
addParameter(p,validatestring('method',{'meth','method'}),...
    fscript.method,@ischar)
addParameter(p,'icer',fscript.icer,@ischar);

% insert into structure
parse(p,callingFunction,varargin{:})

% check for sizeUnit, must be 'mum'
if ~(strcmpi(p.Results.sizeunit,'um') || strcmpi(p.Results.sizeunit,'mum'))
    error('unit for sizes of scatterers must be ''um'' or ''mum''')
end

% check whether to use sensor, wavelengths, or bandPass
newWavelength = ~issens && iswave;
newSensor = ~iswave && issens;
if newWavelength
    fscript.sensor = [];
    fscript.wavelength = p.Results.wavelength(:);
elseif newSensor
    if iscategorical(p.Results.sensor)
        fscript.sensor = p.Results.sensor;
    else
        fscript.sensor = categorical({p.Results.sensor});
    end
end

fscript.waveUnit = p.Results.waveunit;
assert(~isempty(fscript.waveUnit),'''waveUnit'' must be specified')

% if sensor specified, convert to wavelengths or bandPass
if ~isempty(fscript.sensor)
    Tbl = SensorTable(fscript.sensor,fscript.waveUnit);
    fscript.sensor = unique(Tbl.Sensor);
    % spectrometer
    if contains(callingFunction,'spectral','IgnoreCase',true) ||...
            contains(char(fscript.sensor),'aviris','IgnoreCase',true) ||...
            contains(char(fscript.sensor),'casi','IgnoreCase',true)
        fscript.wavelength = Tbl.CentralWavelength;
        fscript.spectrometer = true;
        fscript.multispectral = false;
        if contains(Tbl.Properties.VariableNames,'lowerwavelength','IgnoreCase',true)
            fscript.bandPass = [Tbl.LowerWavelength Tbl.UpperWavelength];
        else
            fscript.bandPass = [];
        end
        %multispectral sensor
    elseif contains(callingFunction,'intg','IgnoreCase',true)
        fscript.spectrometer = false;
        fscript.multispectral = true;
        if isempty(p.Results.bands)
            bandPass = [Tbl.LowerWavelength Tbl.UpperWavelength];
            bands = Tbl.Band;
        else
            bTbl = table(categorical(p.Results.bands(:)),'VariableNames',{'Band'});
            nTbl = innerjoin(Tbl,bTbl);
            bandPass = [nTbl.LowerWavelength nTbl.UpperWavelength];
            bands = nTbl.Band;
        end
        fscript.bandPass = bandPass;
        fscript.bands = bands;
    end
    
elseif ~isempty(p.Results.wavelength)
    fscript.wavelength = p.Results.wavelength;
    assert(isvector(fscript.wavelength),...
        '''wavelength'' must be a vector, not a matrix');
    fscript.wavelength = fscript.wavelength(:);
    fscript.spectrometer = true;
    fscript.multispectral = false;
    
elseif isempty(p.Results.bandpass)
    error('either ''wavelength'', ''sensor'', or ''bandPass'' must be specified')
    
else
    fscript.bandPass = p.Results.bandpass;
    fscript.spectrometer = false;
    fscript.multispectral = true;
    if isempty(p.Results.bands)
        fscript.bands = 1:size(fscript.bandPass,1);
    else
        fscript.bands = p.Results.bands;
    end
end

% if spectrum, radius and cosZ and wavelength must be same size
if isfield(fscript,'spectrometer') && fscript.spectrometer
    % conversion actually not necessary, but in place for future mods
    if isempty(p.Results.cosZ)
        fscript.cosZ = [];
        fscript.muS = [];
        if isSSA
            [radius,fscript.viewF,fscript.wavelength] =...
                checkSizes(SSA2radius(p.Results.ssa,S.unitsSize),...
                p.Results.viewF,fscript.wavelength);
        else
            [radius,fscript.viewF,fscript.wavelength] =...
                checkSizes(convertLengthUnits(p.Results.radius,...
                p.Results.sizeunit,S.unitsSize),p.Results.viewF,fscript.wavelength);
        end
    else
        if isSSA
            if isempty(p.Results.muS)
                [fscript.cosZ,fscript.muS,radius,fscript.viewF,fscript.wavelength] =...
                    checkSizes(p.Results.cosZ,p.Results.cosZ,p.Results.viewF,...
                    SSA2radius(p.Results.ssa,S.unitsSize),...
                    fscript.wavelength);
            else
                [fscript.cosZ,fscript.muS,radius,fscript.viewF,fscript.wavelength] =...
                    checkSizes(p.Results.cosZ,p.Results.muS,p.Results.viewF,...
                    SSA2radius(p.Results.ssa,S.unitsSize),...
                    fscript.wavelength);
            end
        else
            if isempty(p.Results.muS)
%                 [fscript.cosZ,fscript.muS,radius,fscript.viewF,fscript.wavelength] =...
%                     checkSizes(p.Results.cosZ,p.Results.cosZ,p.Results.viewF,...
%                     convertLengthUnits(p.Results.radius,p.Results.sizeunit,S.unitsSize),...
%                     fscript.wavelength);
                [fscript.cosZ,fscript.muS,radius,fscript.viewF,fscript.wavelength] =...
                    checkSizes(p.Results.cosZ,p.Results.cosZ,...
                    convertLengthUnits(p.Results.radius,p.Results.sizeunit,S.unitsSize),...
                    p.Results.viewF,fscript.wavelength);
            else
                [fscript.cosZ,fscript.muS,fscript.viewF,radius,fscript.wavelength] =...
                    checkSizes(p.Results.cosZ,p.Results.muS,p.Results.viewF,...
                    convertLengthUnits(p.Results.radius,p.Results.sizeunit,S.unitsSize),...
                    fscript.wavelength);
            end
            
        end
    end
    %if multispectral, radius and cosZ must be the same
else
    if isempty(p.Results.cosZ)
        fscript.cosZ = [];
        fscript.muS = [];
        fscript.viewF = p.Results.viewF;
        if isSSA
            radius = SSA2radius(p.Results.ssa,S.unitsSize);
        else
            radius = p.Results.radius;
        end
    else
        if isSSA
            if isempty(p.Results.muS)
                [fscript.cosZ,fscript.muS,fscript.viewF,radius] = checkSizes(p.Results.cosZ,...
                    p.Results.cosZ,p.Results.viewF,SSA2radius(p.Results.ssa,S.unitsSize));
            else
                [fscript.cosZ,fscript.muS,fscript.viewF,radius] = checkSizes(p.Results.cosZ,...
                    p.Results.muS,p.Results.viewF,SSA2radius(p.Results.ssa,S.unitsSize));
            end
        else
            if isempty(p.Results.muS)
                [fscript.cosZ,fscript.muS,fscript.viewF,radius] = checkSizes(p.Results.cosZ,...
                    p.Results.cosZ,p.Results.viewF,...
                    convertLengthUnits(p.Results.radius,p.Results.sizeunit,S.unitsSize));
            else
                [fscript.cosZ,fscript.muS,fscript.viewF,radius] = checkSizes(p.Results.cosZ,...
                    p.Results.muS,p.Results.viewF,...
                    convertLengthUnits(p.Results.radius,p.Results.sizeunit,S.unitsSize));
            end
        end
    end
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
    assert(all(R0.wavelength>0) && all(R0.reflectance(:)>=0) &&...
        all(R0.reflectance(:)<=1),'R0 table bad wavelength or reflectance')
end
fscript.R0 = R0;
fscript.sizeUnit = S.unitsSize;
if isempty(p.Results.lap)
    fscript.includeLAP = false;
else
    fscript.includeLAP = true;
    if ischar(p.Results.lap)
        fscript.LAP = {p.Results.lap};
        nLAP = 1;
    elseif istable(p.Results.lap)
        fscript.LAP = p.Results.lap;
        fscript.lookupMie = false;
        nLAP = 1;
        % check wavelength units
        assert(~isempty(fscript.LAP.Properties.VariableUnits),...
            'if LAP properties as table of refractive indices, units for wavelength must be specified')
        fscript.LAP.wavelength = convertLengthUnits(fscript.LAP.wavelength,...
            fscript.LAP.Properties.VariableUnits{1},fscript.waveUnit);
        fscript.LAP.Properties.VariableUnits = {fscript.waveUnit,''};
        % make sure imaginary part of refractive index is positive
        % (although the Mie routines check for this)
        if imag(fscript.LAP.RefractiveIndex(1))<0
            rx = conj(fscript.LAP.RefractiveIndex);
            fscript.LAP.RefractiveIndex = rx;
        end
    else
        assert(iscell(p.Results.lap),'if LAP is neither char nor table, must be cell')
        fscript.LAP = p.Results.lap;
        nLAP = length(fscript.LAP);
    end
    % LAP fraction, same as # of LAPs
    fscript.LAPfraction = p.Results.lfraction;
    fscript.LAPradius = p.Results.lradius;
    assert(length(fscript.LAPfraction)==nLAP,'size of ''lfraction'' vector must = number of LAPs')
    % if cell vector, specify size values if not default
    if iscell(fscript.LAP)
        if isempty(fscript.LAPradius)
            fscript.LAPradius = nan(1,nLAP);
        end
        for k=1:nLAP
            if isnan(fscript.LAPradius(k))
                switch fscript.LAP{k}
                    case 'dust'
                        fscript.LAPradius(k) = S.defaultDustRadius;
                    case 'soot'
                        fscript.LAPradius(k) = S.defaultSootRadius;
                end
            else
                switch fscript.LAP{k}
                    case 'dust'
                        fscript.LAPradius(k) = resetLimits('dustRadius',...
                            fscript.LAPradius(k),S.dustRadius);
                    case 'soot'
                        fscript.LAPradius(k) = resetLimits('sootRadius',...
                            fscript.LAPradius(k),S.sootRadius);
                end
            end
        end
    elseif istable(fscript.LAP)
        if ~isfield(fscript,'LAPradius')
            fscript.LAPradius = S.defaultDustRadius;
        end
    end
end

% snow
if fscript.substance==snow
    if isempty(p.Results.wetness)
        fscript.wetSnow = false;
    else
        fscript.wetSnow = p.Results.wetness>0;
        if ~isempty(p.Results.wetness)
            fscript.wetness = p.Results.wetness;
        end
    end
    if isinf(p.Results.waterequivalent)
        fscript.WE = p.Results.waterequivalent;
        fscript.deepSnow = true;
    else
        fscript.deepSnow = false;
        fscript.WE = p.Results.waterequivalent;
    end
elseif contains(char(fscript.substance),'cloud','IgnoreCase',true)
    assert(~isempty(p.Results.waterequivalent) &&...
        ~isinf(p.Results.waterequivalent),...
        'if substance is some sort of cloud, you must specify a non-infinite ''waterEquivalent'' in mm')
    fscript.WE = p.Results.waterequivalent;
    if fscript.substance==mixedCloud
        assert(~isempty(p.Results.wetness),...
            'for ''mixedCloud'' you must also specify ''wetness''')
        fscript.wetness = p.Results.wetness;
    end
else
    error(['substance=' fscript.substance ' not recognized']); % shouldn't reach
end

% check that R0 and wavelengths match
if istable(fscript.R0)
    % make sure wavelength units for R0 match those for wavelength
    assert(~isempty(fscript.R0.Properties.VariableUnits{1}),...
        'if R0 is a table, wavelength units must be specified')
    fscript.R0.wavelength = convertLengthUnits(fscript.R0.wavelength,...
        fscript.R0.Properties.VariableUnits{1},fscript.waveUnit);
    fscript.R0.Properties.VariableUnits{1} = fscript.waveUnit;
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

%Mie lookup
fscript.lookupMie = logical(p.Results.lookup);

%refractive index
fscript.icer = p.Results.icer;

%elevation
fscript.elevation = p.Results.elevation;

% add SSA or radius to the output and check sizes
if isfield(fscript,'iceRadius')
    if ~isempty(fscript.iceRadius)
        fscript.SSA = radius2SSA(fscript.iceRadius,fscript.sizeUnit);
    end
elseif isfield(fscript,'SSA')
    if ~isempty(fscript.SSA)
        fscript.iceRadius = SSA2radius(fscript.SSA,fscript.sizeUnit);
    end
end

% check sizes against limits
switch fscript.substance
    case snow
        fscript.iceRadius = resetLimits('snowRadius',...
            fscript.iceRadius,S.snowRadius);
        fscript.SSA = resetLimits('snowSSA',fscript.SSA,S.snowSSA);
    case iceCloud
        fscript.iceRadius = resetLimits('iceCloudRadius',...
            fscript.iceRadius,S.iceCloudRadius);
        fscript.SSA = resetLimits('iceCloudSSA',fscript.SSA,S.iceCloudSSA);
    case waterCloud
        fscript.waterRadius = resetLimits('waterCloudRadius',...
            fscript.waterRadius,S.waterCloudRadius);
    case mixedCloud
        fscript.iceRadius = resetLimits('iceCloudRadius',...
            fscript.iceRadius,S.iceCloudRadius);
        fscript.SSA = resetLimits('iceCloudSSA',fscript.SSA,S.iceCloudSSA);
        fscript.waterRadius = resetLimits('waterCloudRadius',...
            fscript.waterRadius,S.waterCloudRadius);
end

% warn for R0 when WE is infinite
if isinf(fscript.WE)
    if isobject(fscript.R0)
        warning('''R0'' is specified even though water equivalent is Inf, so ''R0'' is ignored')
    elseif isscalar(fscript.R0) && fscript.R0>0
        warning('''R0'' is specified even though water equivalent is Inf, so ''R0'' is ignored')
    end
end
end

function outX=resetLimits(iden,inX,XLimits)
t = inX<min(XLimits)*(1-1000*eps);
if any(t)
    warning('variable %s (%g) < minimum (%g), reset',...
        iden,min(inX),min(XLimits))
    inX(t) = min(XLimits)*(1-1000*eps);
end
t = inX>max(XLimits)*(1+1000*eps);
if any(t)
    warning('variable %s (%g) > maximum (%g), reset',...
        iden,max(inX),max(XLimits))
    inX(t) = max(XLimits)*(1+1000*eps);
end
outX = inX;
end