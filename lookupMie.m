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
%Output
%       omega - single scattering albedo
%       g - Mie asymmetry parameter%   M - structure with variables
%       Qext - extinction efficiency
%       Qabs - absorption efficiency
%       Qsca - scattering efficiency
%       Qpr - radiation pressure efficiency

persistent Fice Fwater Fdust Fsoot FwetSnow mievariable

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
            mievariable = S.variable;
        end
        F = Fice;
    case 'water'
        if isempty(Fwater)
            S = load(filename);
            Fwater = S.F;
            mievariable = S.variable;
        end
        F = Fwater;
    case 'dust'
        if isempty(Fdust)
            S = load(filename);
            Fdust = S.F;
            mievariable = S.variable;
        end
        F = Fdust;
    case 'soot'
        if isempty(Fsoot)
            S = load(filename);
            Fsoot = S.F;
            mievariable = S.variable;
        end
        F = Fsoot;
    case 'wetsnow'
        wetSnow = true;
        if isempty(FwetSnow)
            S = load(filename);
            FwetSnow = S.F;
            mievariable = S.variable;
        end
        F = FwetSnow;
    otherwise
        error('''substance'' %s unrecognized',p.Results.substance)
end

if wetSnow
    [radius,lambda,waterFraction] =...
        checkSizes(p.Results.radius,p.Results.lambda,p.Results.waterFraction);
else
    [radius,lambda] = checkSizes(p.Results.radius,p.Results.lambda);
end

%lookup tables use 'mum' for wavelength and radius
sizeUnit = 'um';
waveUnit = 'um';
radius = convertUnits(radius,p.Results.radiusUnits,sizeUnit);
lambda = convertUnits(lambda,p.Results.lambdaUnits,waveUnit);

% raw lookup functions use sqrt of radius and log of wavelength
radius = sqrt(radius);
lambda = log(lambda);

% populate output structure using interpolation function
for k=1:length(mievariable)
    thisF = F{k};
    % for omega near 1.0, lookup uses 1-omega
    if strcmp(mievariable{k},'omega') &&...
            (contains(p.Results.substance,'ice','IgnoreCase',true) ||...
            contains(p.Results.substance,'snow','IgnoreCase',true) ||...
            contains(p.Results.substance,'water','IgnoreCase',true))
        if wetSnow
            M.(mievariable{k}) = 1-thisF(waterFraction,radius,lambda);
        else
            M.(mievariable{k}) = 1-thisF(radius,lambda);
        end
    else
        if wetSnow
            M.(mievariable{k}) = thisF(waterFraction,radius,lambda);
        else
            M.(mievariable{k}) = thisF(radius,lambda);
        end
    end
end

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