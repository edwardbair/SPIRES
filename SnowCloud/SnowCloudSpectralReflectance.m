function [refl,varargout ] = SnowCloudSpectralReflectance(cosZ,substance,radius, radiusUnits, lambda, lambdaUnits, varargin)
% [refl ] = SnowCloudSpectralReflectance(cosZ,substance,radius,radiusUnits,lambda,lambdaUnits [, ...)
% [refl,wave ] = SnowCloudSpectralReflectance(cosZ,substance,radius,radiusUnits,lambda,lambdaUnits [, ...)
% [refl,wave,trans ] = SnowCloudSpectralReflectance(cosZ,substance,radius,radiusUnits,lambda,lambdaUnits [, ...)
%spectral reflectance, and optionally transmittance, of snow
%
%Inputs (if not scalars, must be same size)
% cosZ - cosine of illumination angle
% substance - 'ice' or 'snow' (same), 'water', 'wetSnow', or
% 'mixedPhaseCloud'
%   (to model a mixed phase cloud, specifiy 'mixedPhaseCloud' as the
%   substance, iceRadius as the radius, 'waterConc' as the mass fraction of
%   water, and 'cloudRadius' as the size of the water droplets)
% radius - either scalar or vector of size(lambda)
% radiusUnits - any metric length ('angstrom', 'nm', 'um' or 'mum', 'mm', 'cm', 'm')
% lambda - wavelength of light
% lambdaUnits - any length like radiusUnits, or 'MHz', 'GHz', 'THz'
%
%Optional input, name-value pairs in any order
% 'WE' - snow or cloud water equivalent (Inf if not specified)
% 'waterConc' - mass fraction of water, needed  if 'wetSnow' or
%   'mixedPhaseCloud' is the substance
% 'weUnits' - any metric length (typically 'mm', 'cm', or 'm'), default 'mm'
% 'R0' - reflectance of underlying ground (ignored unless WE specified),
%   can be scalar or same size as lambda, default 0.0,
% 'dust' - dust concentration (mass fraction)
% 'dustRadius' = optical radius of dust
% 'soot' - soot concentration (mass fraction)
% 'sootRadius' - optical radius of soot
% 'cloudRadius' - needed if substance is 'mixedPhaseCloud', same units as
%   radius
% 'lookup', if true, lookup tables are called for Mie calculations if they
%   exist, with a warning if they don't or if input values are out of range,
%   if false, Mie calculations are done (default true) (
%
%Output
% refl - snow or cloud spectral reflectance, dimensionless, same size
%   dimensions as radius or lambda
%Optional output
% wave - wavelengths for which values are calculated (might be different than
%   input lambda in those wavelengths where the current Mie routines do not
%   return correct results (i.e. where real part of refractive index is
%   less than 1.0)
% trans - snow or cloud spectral transmittance (both direct and diffuse)

narginchk(6,22)
nargoutchk(0,3)

warning('%s is an older, deprecated version of my code and will be removed in the future. Use SnowCloudSpectralRefl in my toolbox\SnowCloudReflectance instead.',mfilename)

% parse inputs (converts radii to mum, wavelengths to nm, WE to mm
optargin=size(varargin,2);
assert (mod(optargin,2)==0,'must be even number of optional arguments')
iStruct = parseInput(cosZ,substance,radius,radiusUnits,lambda,lambdaUnits,varargin{:});

% Mie calculations
if iStruct.mixedPhase
    if iStruct.lookup % Mie parameters by lookup table
        M = lookupMie('ice',iStruct.radius,iStruct.sizeUnit,...
            iStruct.lambda,iStruct.waveUnit);
        D(1,1) = lookupMie('water',iStruct.cloudRadius,iStruct.sizeUnit,...
            iStruct.lambda,iStruct.waveUnit);
    else
        M = MieSphere(iStruct.radius,iStruct.sizeUnit,iStruct.lambda,...
            iStruct.waveUnit,'substance','ice');
        D(1,1) = MieSphere(iStruct.cloudRadius,iStruct.sizeUnit,iStruct.lambda,...
            iStruct.waveUnit,'substance','water');
    end
    density = [917 1000];
    AllRad = [iStruct.radius iStruct.cloudRadius];
    contamConc(:,1) = iStruct.waterConc;
    k = 2;
elseif iStruct.wetSnow
    if iStruct.lookup
        M = lookupMie(iStruct.substance,iStruct.radius,iStruct.sizeUnit,...
            iStruct.lambda,iStruct.waveUnit,iStruct.waterConc);
    else
        CRefIn = RefractiveIndex(iStruct.lambda,'ice',iStruct.waveUnit).*(1-iStruct.waterConc) +...
            RefractiveIndex(iStruct.lambda,'water',iStruct.waveUnit).*iStruct.waterConc;
        M = MieSphere(iStruct.radius,iStruct.sizeUnit,iStruct.lambda,...
            iStruct.waveUnit,'refindex',CRefIn);
    end
    density = 917*(1-iStruct.waterConc)+1000*iStruct.waterConc;
    AllRad = iStruct.radius;
    k = 1;
else
    if iStruct.lookup
        M = lookupMie(iStruct.substance,iStruct.radius,iStruct.sizeUnit,...
            iStruct.lambda,iStruct.waveUnit);
    else
        M = MieSphere(iStruct.radius,iStruct.sizeUnit,iStruct.lambda,...
            iStruct.waveUnit,'substance',iStruct.substance);
    end
    if contains(iStruct.substance,'ice','IgnoreCase',true) ||...
            contains(iStruct.substance,'snow','IgnoreCase',true)
        density = 917;
    else
        density = 1000;
    end
    k = 1;
    AllRad = iStruct.radius;
end
if iStruct.dustySnow
    if iStruct.lookup
        D(k) = lookupMie('dust',iStruct.dustRadius,iStruct.sizeUnit,...
            iStruct.lambda,iStruct.waveUnit);
    else
        D(k) = MieSphere(iStruct.dustRadius,iStruct.sizeUnit,iStruct.lambda,...
            iStruct.waveUnit,'substance','dust');
    end
    density = cat(2,density,repmat(2800,size(density,1),1));
    contamConc(:,k) = iStruct.dustConc;
    AllRad = cat(2,AllRad,iStruct.dustRadius);
    k = k+1;
end
if iStruct.sootySnow
    if iStruct.lookup
        D(k) = lookupMie('soot',iStruct.sootRadius,iStruct.sizeUnit,...
            iStruct.lambda,iStruct.waveUnit);
    else
        D(k) = MieSphere(iStruct.sootRadius,iStruct.sizeUnit,iStruct.lambda,...
            iStruct.waveUnit,'substance','soot');
    end
    density = cat(2,density,repmat(2800,size(density,1),1));
    contamConc(:,k) = iStruct.sootConc;
    AllRad = cat(2,AllRad,iStruct.sootRadius);
    k = k+1;
end

% averaged Mie parameters if dirty snow, or mixed phase
% cloud
if ~iStruct.cleanSnowCloud || iStruct.mixedPhase
    allM = cell(1,length(M)+length(D));
    allM{1} = M;
    n = 2;
    for k=1:length(D)
        allM{n} = D(k);
        n = n+1;
    end
    conc = [1-sum(contamConc,2) contamConc];
    mix = MieMixture(allM,AllRad,density,conc);
    M = mix;
end

% eliminate some because of incorrect Mie values
badMie = M.omega>1 | M.omega<0 | M.g>1 | M.g<0;
if nnz(badMie)
    fn = fieldnames(M);
    for k=1:length(fn)
        M.(fn{k})(badMie) = [];
    end
    fn = fieldnames(iStruct);
    for k=1:length(fn)
        if ~isscalar(iStruct.(fn{k})) &&...
                length(iStruct.(fn{k}))==length(badMie)
            iStruct.(fn{k})(badMie) = [];
        end
    end
end

% vector/matrix sizes of mu0 and omega must be equal
[mu0,omega] = checkSizes(iStruct.mu0,M.omega);

if all(isinf(iStruct.WE))
    tau = Inf(size(omega));
else
    switch lower(iStruct.substance)
        case {'mixedphase','mixedphasecloud'}
            mfrac = [1-iStruct.waterConc iStruct.waterConc];
            vfrac = mfrac./[917 1000];
            vfrac = vfrac/sum(vfrac);
            [ri,rw,WE,Qext,wQext] =...
                checkSizes(convertUnits(iStruct.radius,iStruct.sizeUnit,'m'),...
                convertUnits(iStruct.cloudRadius(k),iStruct.sizeUnit,'m'),...
                iStruct.WE,M.Qext,D(1).Qext);
            tau = tauSnow(ri,WE,Qext)*vfrac(1)+tauCloud(rw,WE,wQext)*vfrac(2);
        case {'ice','snow','wetsnow'}
            [r,WE,Qext] = checkSizes(convertUnits(iStruct.radius,iStruct.sizeUnit,'m'),iStruct.WE,M.Qext);
            tau = tauSnow(r,WE,Qext);
        case 'water'
            [r,WE,Qext] =...
                checkSizes(convertUnits(iStruct.radius,iStruct.sizeUnit,'m'),...
                iStruct.WE,M.Qext);
            tau = tauCloud(r,WE,Qext);
        otherwise
            error('''substance'' ''%s'' not recognized',iStruct.substance)
    end
end

[g,R0,~] = checkSizes(M.g,iStruct.R0,omega);
[refl,trans] = twostream(mu0,omega,g,tau,'R0',R0);
if iStruct.changeSize
    refl = reshape(refl,iStruct.origSize);
    if nargout>1
        trans = reshape(trans,iStruct.origSize);
    end
end
if nargout>1
    varargout{1} = convertUnits(iStruct.lambda,iStruct.waveUnit,iStruct.origWaveUnit);
    if nargout>2
    varargout{2} = trans;
    end
end

end

function iStruct = parseInput(cosZ,substance,radius,radiusUnits,lambda,lambdaUnits,varargin)
%parse input values to produce input/output raster references

% some parameters
sizeUnit = 'mum';
waveUnit = 'nm';
S = SnowCloudLimits;
defaultWE = Inf;
defaultWEunits = 'mm';
defaultR0 = 0;

p = inputParser;
rangeValidation = @(x) isnumeric(x) && all(x(:)>=0 & x(:)<=1);
positiveValidation = @(x) isnumeric(x) && all(x(:)>0);
addRequired(p,'cosZ',rangeValidation)
addRequired(p,'substance',@ischar);
addRequired(p,'radius',@isnumeric)
addRequired(p,'radiusUnits',@ischar)
addRequired(p,'lambda',positiveValidation)
addRequired(p,'lambdaUnits',@ischar)
addParameter(p,'we',defaultWE,positiveValidation)
addParameter(p,'weunits',defaultWEunits,@ischar)
addParameter(p,'r0',defaultR0,rangeValidation)
addParameter(p,'dust',[],rangeValidation)
addParameter(p,'dustradius',S.defaultDustRadius,positiveValidation)
addParameter(p,'soot',[],rangeValidation)
addParameter(p,'sootradius',S.defaultSootRadius,positiveValidation)
addParameter(p,'cloudradius',S.defaultWaterCloudRadius,positiveValidation)
addParameter(p,'lookup',true,@islogical)
addParameter(p,'waterconc',[],rangeValidation)
parse(p,cosZ,substance,radius,radiusUnits,lambda,lambdaUnits,varargin{:})

% some logical variables
iStruct.wetSnow = contains(p.Results.substance,'wetsnow','IgnoreCase',true);
iStruct.mixedPhase = (contains(p.Results.substance,'mixedphase','IgnoreCase',true));
iStruct.dustySnow = ~isempty(p.Results.dust);
iStruct.sootySnow = ~isempty(p.Results.soot);
iStruct.deepSnow = isinf(p.Results.we);
iStruct.cleanSnowCloud = ~(iStruct.dustySnow || iStruct.sootySnow);

% units
iStruct.sizeUnit = sizeUnit;
iStruct.waveUnit = waveUnit;
iStruct.origWaveUnit = p.Results.lambdaUnits;

% begin transfer to data structure
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
[mu0,radius,R0,lambda] = checkSizes(p.Results.cosZ,p.Results.radius,...
    p.Results.r0,p.Results.lambda);
if iscolumn(lambda)
    iStruct.changeSize = false;
else
    iStruct.changeSize = true;
    iStruct.origSize = size(lambda);
end
iStruct.mu0 = mu0(:);
iStruct.radius = convertUnits(radius(:),p.Results.radiusUnits,sizeUnit);
iStruct.R0 = R0(:);
iStruct.lambda = convertUnits(lambda(:),p.Results.lambdaUnits,waveUnit);

if iStruct.deepSnow
    iStruct.WE = p.Results.we;
else
    iStruct.WE = convertUnits(p.Results.we,p.Results.weunits,'mm');
end

if ~iStruct.cleanSnowCloud
    if iStruct.dustySnow
        dustConc = p.Results.dust;
        if p.Results.dustradius == S.defaultDustRadius
            dustRadius = convertUnits(p.Results.dustradius,S.unitsSize,iStruct.sizeUnit);
        else
            dustRadius = convertUnits(p.Results.dustradius,...
                p.Results.radiusUnits,iStruct.sizeUnit);
        end
        [iStruct.dustConc,iStruct.dustRadius,~] =...
            checkSizes(dustConc,dustRadius,radius);
    end
    if iStruct.sootySnow
        sootConc = p.Results.soot;
        if p.Results.sootradius == S.defaultSootRadius
            sootRadius = convertUnits(p.Results.sootradius,S.unitsSize,iStruct.sizeUnit);
        else
            sootRadius = convertUnits(p.Results.sootradius,...
                p.Results.radiusUnits,iStruct.sizeUnit);
        end
        [iStruct.sootConc,iStruct.sootRadius,~] =...
            checkSizes(sootConc,sootRadius,radius);
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
                p.Results.radiusUnits,iStruct.sizeUnit);
        end
        [iStruct.waterConc,iStruct.cloudRadius,~] = checkSizes(waterConc,cloudRadius,radius);
    else
        [iStruct.waterConc,~] = checkSizes(waterConc,radius);
    end
end

% use Mie Lookup tables
iStruct.lookup = p.Results.lookup;
end