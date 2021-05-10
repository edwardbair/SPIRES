function [oStruct,varargout] = solveSnowCloudProperties(cosZ,substance,sensor,reflectance,varargin )
% [oStruct] = solveSnowCloudProperties(cosZ,substance,sensor,reflectance,... )
% [oStruct,stats] = solveSnowCloudProperties(cosZ,substance,sensor,reflectance,... )
%
%solves for snow or cloud properties based on spectrum
%properties to solve for depend on substance
%   'snow' - snow, solve for optical grain radius (can specify 'WE'
%       for shallow snow, 'waterFraction' for wet snow, or 'dust' or 'soot'
%       properties)
%   'wetSnow' - grain radius and water fraction (can specify 'WE' for shallow
%       snow, 'dust' concentration for dusty snow, 'soot' concentration for
%       sooty snow)
%   'dustySnow' - solve for grain radius, dust concentration, dust radius (can also be
%       specified)
%   'sootySnow' - solve for grain radius, soot concentration, soot radius (can also be
%       specified)
%   'iceCloud' - solve for grain radius, water equivalent (mm)
%   'waterCloud' - solve for droplet radius, water equivalent (mm)
%   'mixedPhaseCloud' - solve for droplet radius, ice crystal radius,
%       water fraction, water equivalent (mm)
%
%Input
%   cosZ - cosine of solar illumination angle (scalar)
%   substance - 'snow', 'wetSnow', 'dustySnow', 'sootySnow', 'iceCloud',
%       'waterCloud', or 'mixedPhaseCloud' (character string)
%   sensor - any from SensorTable, run SensorTable('help') for list
%   reflectance - vector of reflectance values at angle arccos(cosZ), of
%       length corresponding to number of bands in sensor
%
%Optional input, name-value pairs, case insensitive
%   'band' - either numeric vector, cell vector, or categorical vector of bands,
%       or, if omitted, all bands for that sensor
%   'sizeUnit' - units for radii, default um
%   'WE' - water equivalent for shallow snow, mm (don't specify any of the
%       cloud options)
%   'dust' - dust concentration (don't specify if 'dustySnow')
%   'dustRadius' - can specify, or solve for
%   'soot' - soot concentration (don't specify if 'sootySnow')
%   'sootRadius' - can specify, or solve for
%   'R0' - reflectance of substrate, or if 'fractional' is specified, reflectance
%       of background either scalar or for all bands (either
%       specified by 'band' or length of sensor table), default 0
%   'method' - solution method, default 'lsqnonlin', options 'fminsearch'
%       or 'fmincon' -- NOT YET IMPLEMENTED, using 'lsqnonlin' only
%   'multipleStart' - followed by number of random starting values (default
%       just 1 set of starts) -- NOT YET IMPLEMENTED
%   'fractional' - logical variable to indicate whether to solve for fSCA,
%       applies only if one of the snow members is set in 'substance'
%
%Output
%   oStruct - snow or cloud properties, depending on inputs
%Optional output
%   stats - statistics about solution

warning('function %s deprecated, use the toolbox/SnowCloudReflectance instead', mfilename)

narginchk(4,26)
nargoutchk(0,2)

% parse inputs (convert radii to mum, wavelengths to nm, WE to mm,
optargin=size(varargin,2);
assert (mod(optargin,2)==0,'must be even number of optional arguments')
iStruct = parseInput(cosZ,substance,sensor,reflectance,varargin{:});

% make sure wavelenghts long enough to get size of scatterer
assert(convertUnits(max(iStruct.bandPass(:)),iStruct.waveUnit,'nm')>=1060,...
    'maximum wavelength must be >= %f %s to retrieve size of snow or cloud scatterer',...
    convertUnits(1060,'nm',iStruct.waveUnit),iStruct.waveUnit)

% set unknowns and bounds
S = SnowCloudLimits;
x0 = zeros(1,length(iStruct.solveFor));
lb = zeros(size(x0));
ub = zeros(size(x0));
if contains(iStruct.substance,'snow','IgnoreCase',true)
    for k=1:length(x0)
        switch lower(iStruct.solveFor{k})
            case 'radius'
                x0(k) = S.defaultSnowRadius;
                lb(k) = S.snowRadius(1);
                ub(k) = S.snowRadius(2);
            case 'waterfraction'
                x0(k) = mean(S.wetSnow);
                lb(k) = S.wetSnow(1);
                ub(k) = S.wetSnow(2);
            case 'dust'
                x0(k) = mean(S.dust);
                lb(k) = S.dust(1);
                ub(k) = S.dust(2);
            case 'dustradius'
                x0(k) = S.defaultDustRadius;
                lb(k) = S.dustRadius(1);
                ub(k) = S.dustRadius(2);
            case 'soot'
                x0(k) = mean(S.soot);
                lb(k) = S.soot(1);
                ub(k) = S.soot(2);
            case 'sootradius'
                x0(k) = S.defaultSootRadius;
                lb(k) = S.sootRadius(1);
                ub(k) = S.sootRadius(2);
            otherwise
                error('solveFor = %s not recognized',iStruct.solveFor{k})
        end
    end
    if iStruct.solveFSCA
        x0(end+1) = 0.5;
        lb(end+1) = 0;
        ub(end+1) = 1;
    end
elseif contains(iStruct.substance,'cloud','IgnoreCase',true)
    for k=1:length(x0)
        switch lower(iStruct.solveFor{k})
            case 'iceradius'
                x0(k) = S.defaultIceCloudRadius;
                lb(k) = S.iceCloudRadius(1);
                ub(k) = S.iceCloudRadius(2);
            case 'waterradius'
                x0(k) = S.defaultWaterCloudRadius;
                lb(k) = S.waterCloudRadius(1);
                ub(k) = S.waterCloudRadius(2);
            case 'we'
                switch lower(iStruct.substance)
                    case 'icecloud'
                        x0(k) = mean(S.iceCloudWE);
                        lb(k) = S.iceCloudWE(1);
                        ub(k) = S.iceCloudWE(2);
                    case {'watercloud','mixedphasecloud'}
                        x0(k) = mean(S.waterCloudWE);
                        lb(k) = S.waterCloudWE(1);
                        ub(k) = S.waterCloudWE(2);
                    otherwise
                        error('substance %s not recognized',iStruct.substance)
                end
            case 'waterfraction'
                x0(k) = .5;
                lb(k) = 0;
                ub(k) = 1;
            otherwise
                error('solveFor = %s not recognized',iStruct.solveFor{k})
        end
    end
else
    error('substance = %s not recognized',iStruct.substance)
end

% fill out the arguments for the knowns
lookFor = {'dust','dustRadius','soot','sootRadius','waterConc'};
argcConst = {'bandPass',iStruct.bandPass,'waveUnit',iStruct.waveUnit};
n = 1+length(argcConst);
fn = fieldnames(iStruct);
for k=1:length(lookFor)
    if any(strcmpi(fn,lookFor{k}))
        argcConst{n} = lookFor{k};
        argcConst{n+1} = iStruct.(lookFor{k});
        n = n+2;
    end
end

% solving method depends on input
if length(x0)==1 || (iStruct.solveFSCA && length(x0)==2) % snow, just the radius, + maybe fSCA
    % one unknown, one measurement, so exact solution
    if numel(iStruct.reflectance)==1 && size(iStruct.bandPass,1)==1 && ~iStruct.solveFSCA
        [x,fval,exitflag,output] = fzero(@easySnow,x0);
        iStruct.radius = x;
        stats.fval = fval;
        stats.exitflag = exitflag;
        stats.output = output;
    elseif numel(iStruct.reflectance)==size(iStruct.bandPass,1) % least-squares solution
        [x,resnorm,residual,exitflag,output,lam,jacobian] =...
            lsqnonlin(@easySnow,x0,lb,ub);
        if iStruct.solveFSCA
            iStruct.radius = x(1);
            iStruct.fSCA = x(2);
        else
            iStruct.radius = x;
        end
        stats.resnorm = resnorm;
        stats.residual = residual;
        stats.exitflag = exitflag;
        stats.output = output;
        stats.lambda = lam;
        stats.jacobian = jacobian;
    end
else
    if contains(iStruct.substance,'snow','IgnoreCase',true) % all snow
        % make sure wavelenghts short enough to get particulates
        if any(strcmpi(iStruct.solveFor,'dust')) ||...
                any(strcmpi(iStruct.solveFor,'soot'))
            assert(convertUnits(min(iStruct.bandPass(:)),lambdaUnits,'nm')<=700,...
                'minimum wavelength must be <= %f %s to retrieve dust or soot concentration',...
                convertUnits(700,'nm',lambdaUnits),lambdaUnits)
        end
        [x,resnorm,residual,exitflag,output,lam,jacobian] =...
            lsqnonlin(@anySnow,x0,lb,ub);
    elseif contains(iStruct.substance,'mixedPhaseCloud','IgnoreCase',true)
        [x,resnorm,residual,exitflag,output,lam,jacobian] =...
            lsqnonlin(@mixedCloud,x0,lb,ub);
    else % must be iceCloud or waterCloud
        [x,resnorm,residual,exitflag,output,lam,jacobian] =...
            lsqnonlin(@simpleCloud,x0,lb,ub);
    end
    for k=1:length(iStruct.solveFor)
        iStruct.(iStruct.solveFor{k}) = x(k);
    end
    if iStruct.solveFSCA
        iStruct.fSCA = x(end);
    end
    stats.resnorm = resnorm;
    stats.residual = residual;
    stats.exitflag = exitflag;
    stats.output = output;
    stats.lambda = lam;
    stats.jacobian = jacobian;
end
if nargout>1
    varargout{1} = stats;
end
oStruct = rmfield(iStruct,'solveFSCA');
if isfield(oStruct,'waterConc')
    if oStruct.waterConc==0
        oStruct = rmfield(oStruct,'waterConc');
    end
end

    function diffR = easySnow(x)
        refl = SnowCloudReflectance(iStruct.cosZ,'snow',x(1),iStruct.sizeUnit,argcConst{:});
        if iStruct.solveFSCA
            refl = generateMixSpectrum([refl(:) iStruct.R0(:)],[x(end) 1-x(end)]);
        end
        diffR = refl(:)-iStruct.reflectance(:);
    end

    function diffR = anySnow(x)
        % other arguments to SnowCloudReflectance
        argc = buildArgc(iStruct,x,argcConst{:});
        refl = SnowCloudReflectance(iStruct.cosZ,'snow',x(1),iStruct.sizeUnit,argc{:});
        if iStruct.solveFSCA
            refl = generateMixSpectrum([refl(:) iStruct.R0(:)],[x(end) 1-x(end)]);
        end
        diffR = refl(:)-iStruct.reflectance(:);
    end

    function diffR = simpleCloud(x)
        % other arguments to SnowCloudReflectance
        argc = buildArgc(iStruct,x,argcConst{:});
        if contains(iStruct.substance,'ice','IgnoreCase',true)
            thisSubstance = 'ice';
        else
            thisSubstance = 'water';
        end
        refl = SnowCloudReflectance(iStruct.cosZ,thisSubstance,x(1),...
            iStruct.sizeUnit,argc{:});
        diffR = refl(:)-iStruct.reflectance(:);
    end

    function diffR = mixedCloud(x)
        % other arguments to SnowCloudReflectance
        argc = buildArgc(iStruct,x,argcConst{:});
        for ka=1:length(argc)
            if strcmpi(argc{ka},'WaterRadius')
                argc{ka} = 'cloudRadius';
            end
            if strcmpi(argc{ka},'waterFraction')
                argc{ka} = 'waterConc';
            end
        end
        refl = SnowCloudReflectance(iStruct.cosZ,'mixedPhaseCloud',x(1),...
            iStruct.sizeUnit,argc{:});
        diffR = refl(:)-iStruct.reflectance(:);
    end
end

function iStruct = parseInput(cosZ,substance,sensor,reflectance,varargin)
%parse input values

% some parameters
S = SnowCloudLimits;
defaultWE = Inf;
defaultWEunits = S.unitsWE;
defaultR0 = 0;
defaultMethod = 'lsqnonlin';
waveUnit = 'nm';
sizeUnit = 'mum';

p = inputParser;
rangeValidation = @(x) isnumeric(x) && all(x(:)>=0 & x(:)<=1);
validZ = @(x) isnumeric(x) && isscalar(x) && x>0 && x<=1;
positiveValidation = @(x) isnumeric(x) && isscalar(x) && x>0;
bandValidation = @(x) (isrow(x) || iscolumn(x)) &&...
    ((isnumeric(x) && all(x(:)>0)) || iscell(x) || iscategorical(x));
addRequired(p,'cosZ',validZ)
addRequired(p,'substance',@ischar)
addRequired(p,'sensor',@ischar)
addRequired(p,'reflectance',rangeValidation)
addParameter(p,'band',[],bandValidation)
addParameter(p,'sizeunit',sizeUnit,@ischar)
addParameter(p,'we',defaultWE,positiveValidation)
addParameter(p,'weunits',defaultWEunits,@ischar)
addParameter(p,'r0',defaultR0,rangeValidation)
addParameter(p,'dust',[],rangeValidation)
addParameter(p,'dustradius',[],positiveValidation)
addParameter(p,'soot',[],rangeValidation)
addParameter(p,'sootradius',[],positiveValidation)
addParameter(p,'lookup',true,@islogical)
addParameter(p,'waterconc',0,rangeValidation)
addParameter(p,'fractional',false,@islogical)
addParameter(p,'multiplestart',1,positiveValidation)
addParameter(p,'method',defaultMethod,@ischar)
parse(p,cosZ,substance,sensor,reflectance,varargin{:})

% check for options not yet implemented
assert(p.Results.multiplestart==1,'''multipleStart'' option not yet implemented')
assert(strcmpi(p.Results.method,defaultMethod),'''method'' option not yet implemented')

% units
iStruct.sizeUnit = sizeUnit;
iStruct.waveUnit = waveUnit;

% other required variables
iStruct.substance = p.Results.substance;
iStruct.cosZ = p.Results.cosZ;
iStruct.reflectance = p.Results.reflectance;
iStruct.WE = convertUnits(p.Results.we,p.Results.weunits,defaultWEunits);

% unknowns to solve for
iStruct.solveFSCA = p.Results.fractional &&...
    contains(p.Results.substance,'snow','IgnoreCase',true);
switch lower(iStruct.substance)
    case 'snow'
        solveFor = {'radius'};
    case 'wetsnow'
        solveFor = {'radius','waterFraction'};
    case 'dustysnow'
        solveFor = {'radius','dust'};
        if isempty(p.Results.dustradius)
            solveFor = [solveFor 'dustRadius'];
        end
    case 'sootysnow'
        solveFor = {'radius','soot'};
        if isempty(p.Results.sootradius)
            solveFor = [solveFor 'sootRadius'];
        end
    case 'icecloud'
        solveFor = {'IceRadius','WE'};
        iStruct.WE = [];
    case 'watercloud'
        solveFor = {'WaterRadius','WE'};
        iStruct.WE = [];
    case 'mixedphasecloud'
        solveFor = {'IceRadius','WaterRadius','WE','waterFraction'};
        iStruct.WE = [];
end
iStruct.solveFor = solveFor;

% begin transfer to data structure
% sensor/band characteristics
iStruct.bandPass = getBands(p,iStruct);

% check and adjust sizes}
[iStruct.R0,~] = checkSizes(p.Results.r0,iStruct.bandPass(:,1));
lookFor = {'dust','dustRadius','soot','sootRadius','waterConc'};
fn = fieldnames(p.Results);
for k=1:length(lookFor)
    if any(strcmpi(fn,lookFor{k}))
        if ~isempty(p.Results.(lower(lookFor{k})))
            iStruct.(lookFor{k}) = p.Results.(lower(lookFor{k}));
        end
    end
end
end

function bandPass = getBands(p,iStruct)
% check consistency of band designations
T = SensorTable(p.Results.sensor,iStruct.waveUnit);
if isempty(p.Results.band)
    % all bands
    if strcmpi(p.Results.sensor,'aviris-ng')
        bandPass = T.CentralWavelength;
    else
        bandPass = [T.LowerWavelength T.UpperWavelength];
    end
else
    x = p.Results.band;
    if isnumeric(x)
        band = categorical(x);
    elseif iscategorical(x)
        band = x;
    else % cell
        band = categorical(x);
    end
    % need to add stuff for spectrometer here
    bandPass = zeros(length(x),2);
    for k=1:length(band)
        b = find(T.Band==band(k));
        if isempty(b)
            warning(['band ' band(k) ' not found'])
        else
            bandPass(k,:) = [T.LowerWavelength(b) T.UpperWavelength(b)];
        end
    end
end
% make sure bigger number is in column 2
if ~isvector(bandPass)
    bandPass = sort(bandPass,2);
end
end

function argc=buildArgc(S,x,varargin)
if nargin>2
    argc = varargin;
    n = length(argc)+1;
else
    n = 1;
end
for k=2:length(S.solveFor) % start w index 2 because radius is #1
    argc{n} = S.solveFor{k};
    argc{n+1} = x(k);
    n = n+2;
end
end