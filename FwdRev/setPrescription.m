function prescription = setPrescription(substance,varargin)
% prescription=setPrescription(substance,prop1,val1,pro2,val2,...)
%setPrescription: defines the prescription for a reflectance model based on
%snow or cloud properties
% A set of property/value pairs are parsed into structures to use as the
% prescription for functions that calculate reflectance, emissivity, and
% transmittance
%
% (inspired by John D'Errico's slmset for the SLM toolbox, available on
% the MATLAB File Exchange)
%
%First argument, required:
% substance,either 'snow', 'iceCloud', 'waterCloud', or 'mixedCloud'
%   (but any unambiguous abbreviation beginning with first letter works)
%
%% remaining arguments
%   prop/val as a set of name-value pairs (detailed descriptions below)
%       property names are case-insensitive and can be shortened as long as
%       the short name is unambiguous
%%  illumination geometry
%   'cosZ' - cosine of illumination angle on surface, scalar, vector, or matrix
%       (default value 2/3, explictly set to empty to calculate diffuse reflectance)
%   'muS' - cosine of illumination angle on slope (default same as cosZ)
%   'viewF' - view factor (default 1)
%   'elevation' - in m
%% Wavelength properties, either 'wavelength', 'bandPass', or 'sensor'/'band'
%   must be specified, but just one
% 'wavelength' – vector if sensor is a spectrometer
% 'bandPass' Nx2 matrix if sensor is multispectral
% 'waveUnit' – units for wavelength, no default (to prevent errors)
% 'sensor' – instead of 'wavelength' or 'bandPass', can specify spectrometer
%   or multispectral sensor, anything in the SensorTable.m function
%   (if sensor is a spectrometer, returns vector of Central Wavelengths, or
%   if multispectral sensor, returns combined vector of Central, Lower, and
%   Upper Wavelengths)
% 'bands' – bands of the sensor, either numeric vector, cell vector, or
%   categorical vector of bands,or, if omitted, all bands for that sensor
%% Properties applicable to snow
% 'sizeUnit' – units for optically equivalent radius of snow grains (default
%   'mum', none others accepted -- might fix sometime)
% 'radius' - radius of optically equivalent sphere
% 'LAP' - light-absorbing particles, followed either by 'dust' or 'soot',
%   or optionally by a table with a column for wavelength (with units for
%   wavelength specified in Tbl.Properties.VariableUnits) and a second
%   column, titled RefractiveIndex, which is a matrix with a column for each
%   LAP, as a complex number with imaginary parts positive or negative
%   (the wavelengths need not be the same as those specified by 'wavelength'
%   or 'sensor', but if they don't cover the range, nearest values will be
%   used for extrapolation)
% 'LAPfraction' – mass fraction of light-absorbing particles, as a vector if
%   multiple LAPS
% 'LAPradius' – in same units as for optically equivalent snow grain radius,
%   as a vector if multiple LAPs (can use default if 'dust' or 'soot' or
%   default is assumed if 'lookup' is true)
% 'waterEquivalent' – water equivalent of snow in mm (same as kg/m2), Inf if
%   not specified
% 'wetness' – water mass fraction (0 to 0.2)
% 'waterRadius' - default 5 um
% 'temperature' - degrees K, physical temperature of snow, used only in
%   calculations involving wavelength-integrated emissivity
% 'fSCA' - fractional snow-covered area (0-1) + a vector of length of the
%   2nd dimension of R0 so that sum(fSCA)=1 and the 2nd through last
%   endmembers have reflectance R0
%% Properties applicable to clouds
% 'sizeUnit' – units for optically equivalent radius of snow grains (default
%   'mum', none others accepted -- might fix sometime)
% 'radius' - radius of optically equivalent sphere
% 'waterEquivalent' – water equivalent of cloud in mm (same as kg/m2)
% 'wetness' – water mass fraction, >0 but <1 if specified, otherwise 0
% 'waterRadius' - default 5 um
% 'temperature' - degrees K, physical temperature of cloud, used only in
%   calculations involving wavelength-integrated emissivity
%% Substrate properties
% 'R0' – reflectance of surface beneath snow or cloud or next to snow if
%   waterEquivalent is Inf (in this case spectral mixing)
%   Can specify as a scalar
%   Can specify as a  griddedInterpolant or cfit object that takes wavelength
%       as the argument (in this case code does not check if its wavelength
%       has the correct units)
%   Can specify as a cell vector of griddedInterpolant or cfit objects
%   Can specify as a table with columns wavelength and reflectance, with
%       units for wavelength specified in Tbl.Properties.VariableUnits
%       (must specify, will convert if different than waveUnit)
%% Properties about the radiative transfer calculations
% 'icer' - which refractive index data, 'icewb' for Warren & Brandt 2008,
%   'icep' for Picard modification to Warren & Brandt 2008, 'icew' for Warren
%   1984 (default is 'icep')
% 'lookup' – true or false to use lookup tables to calculate Mie variables, no default
% 'method' - (default is delta-Eddington, can specify by abbreviations)
%   'meador' or 'hybrid' - for Meader-Weaver hybrid
%   'delta' or 'eddington' - for delta-Eddington (like Wiscombe & Warren)
%   'disort' (discrete ordinates, not yet implemented)
% 'diffuseFraction' - function of wavelength to estimate diffuse fraction
%   (needed for terrain-modified calculations)
% (generally any 3- or more-letter abbreviation works)

%%
narginchk(1,Inf)
nargoutchk(0,1)
S = SnowCloudLimits;
strings = {'snow','ice','iceCloud','water','waterCloud','mixed',...
    'mixedCloud','mixedPhase','mixedPhaseCloud'};

% process inputs (complicated, lots of options)
p = inputParser;
addRequired(p,'substance',@(x) ischar(x) || isstring(x) || iscategorical(x));

%% all substances
snow = categorical({'snow'});
iceCloud = categorical({'iceCloud'});
waterCloud = categorical({'waterCloud'});
mixedCloud = categorical({'mixedCloud'});
notSet = categorical({'notSet'}); %#ok<NASGU>

%% set the substance
matchstr = validatestring(char(substance),strings);
switch matchstr
    case 'snow'
        prescription.substance = snow;
        prescription.snow = struct;
    case {'ice','iceCloud'}
        prescription.substance = iceCloud;
        prescription.iceCloud = struct;
        cloudStruct = 'iceCloud';
    case {'water','waterCloud'}
        prescription.substance = waterCloud;
        prescription.waterCloud = struct;
        cloudStruct = 'waterCloud';
    case {'mixed','mixedCloud','mixedPhase','mixedPhaseCloud'}
        prescription.substance = mixedCloud;
        prescription.mixedCloud = struct;
        cloudStruct = 'mixedCloud';
    otherwise
        error('No match for substance ''%s''',matchstr)
end

%% parse input and check sizes
r0Validation = @(x) (isnumeric(x) && isvector(x) && all(x>=0) && all(x<=1)) ||...
    (isobject(x) && (isa(x,'griddedInterpolant') || isa(x,'cfit'))) ||...
    iscell(x) ||...
    (istable(x) &&...
    ~isempty(any(contains(x.Properties.VariableNames,'wavelength','IgnoreCase',true))) &&...
    ~isempty(any(contains(x.Properties.VariableNames,'reflectance','IgnoreCase',true))));
nonNegValidation = @(x) isnumeric(x) && all(x(:)>=0) && all(x(:)<=1);
positiveValidation = @(x) isnumeric(x) && all(x(:)>0);

%illumination angle(s), elevation
addParameter(p,'cosZ',2/3,nonNegValidation)
addParameter(p,'muS',[],nonNegValidation)
addParameter(p,'viewF',1,nonNegValidation)
addParameter(p,validatestring('elevation',{'elev','elevation'}),3000,...
    @(x) isnumeric(x) && all(x(:)>-500) && all(x(:)<9000))

% wavelength or sensor (can specify just one, not both
addParameter(p,validatestring('wavelength',{'wavel','wavelength'}),...
    [],positiveValidation)
addParameter(p,validatestring('waveunit',{'waveu','waveunit','lambda'}),...
    '',@(x) ischar(x) && ~isempty(x))
addParameter(p,validatestring('sensor',{'sens','sensor'}),...
    '',@(x) ischar(x) || iscategorical(x))
addParameter(p,'bands',[],@(x) ischar(x) || iscategorical(x) || isnumeric(x))
addParameter(p,validatestring('bandpass',{'bandp','bandpass'}),[],...
    @(x) isnumeric(x) && size(x,2)==2)
addParameter(p,'ignoresolar',false,@(x) isscalar(x) &&...
    (islogical(x) || isnumeric(x)))

% snow
if isfield(prescription,'snow')
    addParameter(p,validatestring('sizeunit',{'size','sizeunit'}),...
        'mum',@ischar)
    addParameter(p,'radius',S.defaultSnowRadius,@(x) isnumeric(x) &&...
        all(x(:)>=min(S.snowRadius)) && all(x(:)<=max(S.snowRadius)))
    addParameter(p,'lap',[],@(x) ischar(x) || iscell(x) || istable(x));
    addParameter(p,validatestring('lapfraction',{'lfr','lfrac','lfraction','lapfr','lapfraction'}),...
        [],@(x) isnumeric(x) && all(x>=0) && all(x<1))
    addParameter(p,validatestring('lapradius',{'lr','lrad','lradius','laprad','lapradius'}),...
        [],@(x) isnumeric(x) && all(x>0))
    addParameter(p,validatestring('waterequivalent',{'weq','wequiv','waterequiv','waterequivalent'}),...
        Inf,@(x) isscalar(x) && x>0)
    addParameter(p,validatestring('temperature',{'tem','temp','temperature'}),...
        273.16,@(x) isscalar(x) && x>0)
    addParameter(p,validatestring('wetness',{'wetn','wetness'}),0,...
        @(x) isscalar(x) && x>=min(S.wetSnow) && x<=max(S.wetSnow))
    addParameter(p,'fSCA',[1 0],@(x) isvector(x) && isnumeric(x) &&...
        all(x>=0) && all(x<=1))
elseif isfield(prescription,'iceCloud') || isfield(prescription,'waterCloud') ||...
        isfield(prescription,'mixedCloud')
    addParameter(p,validatestring('sizeunit',{'size','sizeunit'}),...
        'mum',@ischar)
    addParameter(p,validatestring('waterequivalent',{'weq','wequiv','waterequiv','waterequivalent'}),...
        [],@(x) isscalar(x) && x>0 && ~isempty(x))
    addParameter(p,validatestring('temperature',{'tem','temp','temperature'}),...
        273.16,@(x) isscalar(x) && x>0)
    if isfield(prescription,'waterCloud')
        addParameter(p,'radius',S.defaultWaterCloudRadius,@(x) isnumeric(x) &&...
            all(x(:)>=min(S.defaultWaterCloudRadius)) &&...
            all(x(:)<=max(S.defaultWaterCloudRadius)))
    else
        addParameter(p,'radius',S.defaultIceCloudRadius,@(x) isnumeric(x) &&...
            all(x(:)>=min(S.defaultIceCloudRadius)) &&...
            all(x(:)<=max(S.defaultIceCloudRadius)))
    end
    if isfield(prescription,'mixedCloud')
        addParameter(p,validatestring('waterradius',{'water','waterr','waterradius'}),...
            S.defaultWaterCloudRadius,@(x) isscalar(x) && isnumeric(x) &&...
            x>=min(S.defaultWaterCloudRadius) &&...
            x<=max(S.defaultWaterCloudRadius))
        addParameter(p,validatestring('wetness',{'wetn','wetness'}),0,...
            @(x) isscalar(x) && x>=0 && x<1);
    end
end

% radiative-transfer
addParameter(p,validatestring('lookup',{'look','lookup'}),...
    true,@(x) islogical(x) || isnumeric(x))
addParameter(p,validatestring('method',{'meth','method'}),...
    'delta-Eddington',@ischar)
addParameter(p,'icer','icep',@ischar);
addParameter(p,validatestring('diffusefraction',{'diff','difff','diffusefraction'}),...
    [],@(x) strcmpi(class(x),'cfit') ||...
    contains(class(x),'interpolant','IgnoreCase',true))

% substrate
addParameter(p,'r0',[],r0Validation);

%% insert into structure
parse(p,substance,varargin{:})

%% wavelength or sensor structure
prescription.Spectrum = fillSpectrum(p.Results);
%% substrate before snow to see if more than one non-snow endmember
prescription.Substrate = fillSubstrate(prescription.Spectrum,p.Results);
%% fill snow or cloud structure, sizeUnit must be 'mum'
if isfield(prescription,'snow')
    prescription.snow = fillSnow(prescription,S,p.Results);
else
    prescription.(cloudStruct) = fillCloud(prescription.(cloudStruct),...
        p.Results,cloudStruct);
    prescription.(cloudStruct).SizeUnit = p.Results.sizeunit;
end
%% illumination geometry
prescription.Illumination = fillIllumination(prescription,p.Results);
%% radiative transfer
prescription.RadiativeTransfer = fillRT(p.Results);
end
