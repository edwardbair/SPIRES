function [ Tb ] = SnowCloudBandPassBrightnessTemp(physTemp,cosZ,lambda, lambdaUnits, varargin)
% [ Tb ] = SnowCloudBandPassBrightnessTemp(physTemp,cosZ,lambda, lambdaUnits, varargin)
%brightness temperature of snow or cloud integrated over wavelength range
%
%Inputs
% physTemp - physical temperature, Kelvin (scalar)
% cosZ - cosine of viewing angle (scalar)
% lambda - wavelength range of light, vector of length 2
% lambdaUnits - any length like 'um', or 'MHz', 'GHz', 'THz'
%Optional input, name-value pair
% 'radius' scalar
% 'radiusUnits' - any metric length ('angstrom', 'nm', 'um' (default) or 'mum', 'mm', 'cm', 'm')
% 'WE' - snow or cloud water equivalent (Inf if not specified)
% 'weUnits' - any metric length (typically 'mm', 'cm', or 'm')
% 'substance' - 'ice' (default) or 'water'
%
%Output
% Tb - brightness temperature, Kelvin

% parse inputs
defaultWE = Inf;
defaultWEunits = '';
defaultSubstance = 'ice';
defaultRadius = 250;
defaultRadiusUnits = 'um';
% must be scalars (otherwise makes code more obtuse)
validationFcn = @(x) isnumeric(x) && isscalar(x);
p = inputParser;
addRequired(p,'physTemp',validationFcn)
addRequired(p,'cosZ',validationFcn)
addRequired(p,'lambda',@isnumeric)
addRequired(p,'lambdaUnits',@ischar)
addParameter(p,'radius',defaultRadius,validationFcn)
addParameter(p,'radiusunits',defaultRadiusUnits,@ischar)
addParameter(p,'we',defaultWE,@isnumeric)
addParameter(p,'weunits',defaultWEunits,@ischar)
addParameter(p,'substance',defaultSubstance,@ischar)
parse(p,physTemp,cosZ,lambda,lambdaUnits,varargin{:})
lambda = p.Results.lambda;
assert(length(lambda)==2,'lambda must be a vector of length 2')

% convert wavelengths to meters
lambda = convertLengthUnits(lambda,p.Results.lambdaUnits,'m');
passLambda = lambda;

% other variables to pass to functions
passPhysTemp = p.Results.physTemp;
passRunits = p.Results.radiusunits;
passWE = p.Results.we;
passWunits = p.Results.weunits;
passSubstance = p.Results.substance;
passRadius = p.Results.radius;
passCosine = p.Results.cosZ;
radiance = integral(@emittedR,lambda(1),lambda(2));

% convert radiance(s) to brightness temperature(s)
Tmelt = 273.15;
passRadiance = radiance;
Tb = fzero(@radDiff,Tmelt);

% nested functions
    function rad = emittedR(wavelength)
        radBlack = planck(wavelength,passPhysTemp);
        refl = SnowCloudSpectralReflectance(passCosine,passRadius,passRunits,...
            wavelength,'m','WE',passWE,'WEunits',passWunits,...
            'substance',passSubstance);
        rad = radBlack.*(1-refl);
    end

    function radD = radDiff(T)
        Rbright = PlanckIntegral(passLambda(1),passLambda(2),T);
        radD = Rbright-passRadiance;
    end

end