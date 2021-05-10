function [ R ] = bandPassEmittedRadiance(lambda, lambdaUnits, spectralEmissivity, temperature, bandPass)
% [ R ] = bandPassReflectance(lambda, lambdaUnits, spectralReflectance, bandPass, [Name,Value])
%average spectral emitted radiance over wavelength range of band pass
%
%Inputs (lambda and spectralEmissivity must be same size)
% lambda - wavelength range of spectralEmissivity values
% lambdaUnits - any length unit, typically 'um', 'mum', or 'nm'
% spectralEmissivity - vector of same size as lambda
% bandPass - matrix of size Nx2, specifying the wavelength ranges of N band
%   passes, must be in the same units as lambdaUnits
% temperature - of emitting surface, Kelvin
%
%Output
% R - band-averaged radiance values, matrix of size Nx1,
%   units W m^(-2) lambdaUnits^(-1) sr^(-1)

% parse inputs

p = inputParser;
nonNegativeFcn = @(x) isnumeric(x)&&all(x(:)>=0);
bpValidation = @(x) isnumeric(x) && all(x(:)>=0) && size(x,2)==2;
addRequired(p,'lambda',nonNegativeFcn)
addRequired(p,'lambdaUnits',@ischar)
addRequired(p,'spectralEmissivity',nonNegativeFcn)
addRequired(p,'temperature',nonNegativeFcn)
addRequired(p,'bandPass',bpValidation)
parse(p,lambda,lambdaUnits,spectralEmissivity,temperature,bandPass)

% check sizes
lambda = p.Results.lambda;
spectralEmissivity = p.Results.spectralEmissivity;
bandPass = p.Results.bandPass;
assert(iscolumn(lambda) || iscolumn(spectralEmissivity) ||...
    isrow(lambda) || isrow(spectralEmissivity),...
    'lambda and spectralEmissivity must be vectors, not matrices')
assert(isequal(numel(lambda),numel(spectralEmissivity)),...
    'lambda and spectralEmissivity vectors must be same size')
% use column vectors for both
lambda = lambda(:);
spectralEmissivity = spectralEmissivity(:);
passTemp = p.Results.temperature;

% make sure bandPass values are within spectrum
tol = 1.e-6;
if ~all(bandPass(:)>=min(lambda)) || ~all(bandPass(:)<=max(lambda))
    if max(bandPass(:)-min(lambda)) < -tol
        error('check input - some bandPass wavelengths < min(lambda)')
    elseif max(bandPass(:)-max(lambda)) > tol
                error('check input - some bandPass wavelengths > max(lambda)')
    end
end

%convert wavelengths to meters for Planck function
lambda = convertUnits(lambda,p.Results.lambdaUnits,'m');
bandPass = convertUnits(bandPass,p.Results.lambdaUnits,'m');

% interpolation function for spectral emissivity
Femiss = fit(lambda,spectralEmissivity,'pchip');

R = zeros(size(bandPass,1),1);
for k=1:size(bandPass,1)
        R(k) = integral(@integrandR,bandPass(k,1),bandPass(k,2),'ArrayValued',true)/...
            (bandPass(k,2)-bandPass(k,1));
end

% convert back to lambdaUnits^(-1)
R = convertUnits(R,p.Results.lambdaUnits,'m'); % note: arguments inverted because ^(-1)

    function X = integrandR(wavelength)
        emiss = Femiss(wavelength)';
        S = planck(wavelength',passTemp);
        X = (emiss.*S);
    end
end