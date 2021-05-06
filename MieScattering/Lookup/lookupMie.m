function M = lookupMie(substance,radius,radiusUnits,lambda,lambdaUnits,varargin)
% M = lookupMie(substance,radius,radiusUnits,lambda,lambdaUnits [,waterFraction])
%lookupMie look up Mie variables for ice, water, dust, or soot
%
%Input
%   substance - either 'ice' or 'snow' (same), 'water', 'wetSnow',
%       'dust', or 'soot' (various different dusts to be added)
%   radius - optical radius of particle
%   radiusUnits - typically 'um', 'mm', etc
%   lambda - wavelengths
%   lambdaUnits = typically 'nm', 'um'
%   (if both radius and lambda are vectors or matrices, they must be same
%   size)
%
%Optional input needed if 'wetSnow' is specified
%   waterFraction - mass fraction of liquid water
%
%Output - M structure with variables
%       omega - single scattering albedo
%       g - Mie asymmetry parameter
%       Qext - extinction efficiency
% (the following are calculated from omega, g, and Qext)
%       Qabs - absorption efficiency
%       Qsca - scattering efficiency
%       Qpr - radiation pressure efficiency

persistent Fice Fwater Fdust Fsoot FwetSnow varSnow varWater varDust varSoot varWet

p = inputParser;
rangeValidation = @(x) isnumeric(x) && all(x(:)>=0);
positiveValidation = @(x) isnumeric(x) && all(x(:)>0);
addRequired(p,'substance',@ischar)
addRequired(p,'radius',positiveValidation)
addRequired(p,'radiusUnits',@ischar)
addRequired(p,'lambda',positiveValidation)
addRequired(p,'lambdaUnits',@ischar)
addOptional(p,'waterFraction',0,rangeValidation)
parse(p,substance,radius,radiusUnits,lambda,lambdaUnits,varargin{:})

% load lookup table
if strcmpi(p.Results.substance,'snow')
    filename = ['LUT_Mie_' 'ice' '.mat'];
else
    filename = ['LUT_Mie_' lower(p.Results.substance) '.mat'];
end
wetSnow = false;
switch lower(p.Results.substance)
    case {'ice','snow'}
        if isempty(Fice)
            S = load(filename);
            Fice = S.F;
            varSnow = S.outputVariable;
        end
        F = Fice;
        mievariable = varSnow;
    case 'water'
        if isempty(Fwater)
            S = load(filename);
            Fwater = S.F;
            varWater = S.outputVariable;
        end
        F = Fwater;
        mievariable = varWater;
    case 'dust'
        if isempty(Fdust)
            S = load(filename);
            Fdust = S.F;
            varDust = S.outputVariable;
        end
        F = Fdust;
        mievariable = varDust;
    case 'soot'
        if isempty(Fsoot)
            S = load(filename);
            Fsoot = S.F;
            varSoot = S.outputVariable;
        end
        F = Fsoot;
        mievariable = varSoot;
    case 'wetsnow'
        wetSnow = true;
        if isempty(FwetSnow)
            S = load(filename);
            FwetSnow = S.F;
            varWet = S.outputVariable;
        end
        F = FwetSnow;
        mievariable = varWet;
    otherwise
        error('''substance'' %s unrecognized',p.Results.substance)
end

if wetSnow
    [waterFraction,lambda,radius] = checkSizes(p.Results.waterFraction,...
        p.Results.lambda,p.Results.radius);
else
    [lambda,radius] = checkSizes(p.Results.lambda,p.Results.radius);
end

%lookup tables use 'um' for wavelength and radius
sizeUnit = 'um';
waveUnit = 'um';
radius = convertLengthUnits(radius,p.Results.radiusUnits,sizeUnit);
lambda = convertLengthUnits(lambda,p.Results.lambdaUnits,waveUnit);

% raw lookup functions use sqrt of radius and log of wavelength
radius = sqrt(radius);
lambda = log(lambda);

% populate output structure using interpolation function
for k=1:length(mievariable)
    thisF = F{k};
    % lookup uses log(1-omega)
    if strcmp(mievariable{k},'omega')
        if wetSnow
            M.(mievariable{k}) = 1-exp(thisF(waterFraction,lambda,radius));
        else
            M.(mievariable{k}) = 1-exp(thisF(lambda,radius));
        end
    else
        if wetSnow
            M.(mievariable{k}) = thisF(waterFraction,lambda,radius);
        else
            M.(mievariable{k}) = thisF(lambda,radius);
        end
    end
end

%lookup covers just omega, g, and Qext -- calculate the others
M.Qsca = M.Qext.*M.omega;
M.Qabs = M.Qext.*(1-M.omega);
M.Qpr = M.Qext-M.g.*M.Qsca;

%wavelength and radius in original units
M.wavelength = unique(p.Results.lambda(:));
M.waveUnits = p.Results.lambdaUnits;
M.radius = unique(p.Results.radius(:));
M.radiusUnits = p.Results.radiusUnits;

% warn if maybe out of range
for k=1:length(mievariable)
    if any(isnan(M.(mievariable{k})))
        limR = [min(F{k}.GridVectors{1}) max(F{k}.GridVectors{1})].^2;
        limW = exp([min(F{k}.GridVectors{2}) max(F{k}.GridVectors{2})]);
        warning('variable %s (and probably others), %d NaN values, inputs maybe out of range, radius limits (%s) are [%g %g], wavelength limits (%s) are [%g %g], check your inputs',...
            mievariable{k},nnz(isnan(M.(mievariable{k}))),sizeUnit,limR(1),limR(2),waveUnit,limW(1),limW(2));
        break
    end
end

end