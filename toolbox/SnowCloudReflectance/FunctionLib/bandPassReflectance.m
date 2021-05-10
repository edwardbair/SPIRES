function [ R,varargout ] = bandPassReflectance(lambda, lambdaUnit, spectralReflectance, spectralIrradiance, varargin)
% [ R ] = bandPassReflectance(lambda, lambdaUnit, spectralReflectance, spectralIrradiance, 'bandPass', values)
% [ R ] = bandPassReflectance(lambda, lambdaUnit, spectralReflectance, 'sensor', sensorName, 'bands', values)
% [R, bandPass] = bandPassReflectance(___________)
%integration over wavelength ranges from spectral reflectance and irradiance
%
%Inputs
% lambda - wavelengths of spectralReflectance values, can be row or column
%   vector
% lambdaUnit - any length unit, typically 'um', 'mum', or 'nm', of length L
% spectralReflectance - vector of same size as lambda, or matrix of size
%   Lx2 or 2xL where L is the length of lambda and columns/rows 1 and 2 are
%   the direct and diffuse spectral reflectance
% spectralIrradiance - vector of same length (L) as lambda or same length as
%   the optional irradWavelength vector as described below among the optional
%   arguments, or matrix of size Lx2 or 2xL where L is the length of lambda
%   (or the length of irradWavelength) and columns/rows 1 and 2 are
%   the direct and diffuse spectral irradiances (can be empty, in which
%   case solar radiation is ignored and function just returns band-averaged
%   reflectivity)
%Optional input, name-value pair (either 'bandPass' or 'sensor'/'bands'
%   must be specified, but not both
% 'bandPass' - matrix of size Nx2, specifying the wavelength ranges of N
%   band passes, must be in the same units as lambdaUnit
% 'sensor' - followed by character string, any sensor from SensorTable.m
% 'band' - followed by either numeric vector, cell vector, or categorical
%   vector of bands, or, if omitted, all bands for that sensor
%   (numeric vector can be used for those sensors whose bands are designated just by number)
% 'irradWavelength' - wavelengths of irradiance, if different than lambda,
%   but must be same units as lambda and same length as spectralIrradiance
%
%Output
% R - band-averaged reflectance values, column vector of length N
% bandPass (optional) - band passes of the input, as a table

% hang onto previous results for interpolation function of solarIrradiance
persistent Last_irradWavelength Last_spectralIrradiance SolDirF SolDifF

% parse inputs
narginchk(5,9)
nargoutchk(0,2)
p = inputParser;
nonNegativeFcn = @(x) isnumeric(x)&&all(x(:)>=0);
bpValidation = @(x) isnumeric(x) && all(x(:)>=0) && (size(x,2)==2 || size(x,1)==2);
bandValidation = @(x) (isrow(x) || iscolumn(x)) &&...
    ((isnumeric(x) && all(x(:)>0)) || iscell(x) || iscategorical(x));
addRequired(p,'lambda',nonNegativeFcn)
addRequired(p,'lambdaUnit',@ischar)
addRequired(p,'spectralReflectance',nonNegativeFcn)
addRequired(p,'spectralIrradiance',nonNegativeFcn)
addParameter(p,'bandpass',[],bpValidation)
addParameter(p,'sensor','',@ischar)
addParameter(p,'band',[],bandValidation)
addParameter(p,'irradWavelength',[],nonNegativeFcn)
parse(p,lambda,lambdaUnit,spectralReflectance,spectralIrradiance,varargin{:})

% check sizes & make column vectors and Lx2 matrices
assert(isrow(p.Results.lambda) || iscolumn(p.Results.lambda),...
    'lambda must be row or column vector')
lambda = p.Results.lambda(:);
spectralReflectance = p.Results.spectralReflectance;
assert(length(lambda)==size(spectralReflectance,1) ||...
    length(lambda)==size(spectralReflectance,2),...
    'lambda and spectralReflectance must be same length')
if isrow(spectralReflectance)
    spectralReflectance = spectralReflectance(:);
elseif length(lambda)==size(spectralReflectance,2)
    spectralReflectance = spectralReflectance';
end
ignoreSolar = isempty(p.Results.spectralIrradiance);
if ~ignoreSolar
    if isempty(p.Results.irradWavelength)
        irradWavelength = lambda;
        spectralIrradiance = p.Results.spectralIrradiance;
        assert(numel(spectralIrradiance)==numel(spectralReflectance),...
            'size of spectralIrradiance must be same as size of spectralReflectance')
        if isrow(spectralIrradiance)
            spectralIrradiance = spectralIrradiance(:);
        elseif length(lambda)==size(spectralIrradiance,2)
            spectralIrradiance = spectralIrradiance';
        end
    else
        irradWavelength = p.Results.irradWavelength;
        assert(iscolumn(irradWavelength) || isrow(irradWavelength),...
            'irradWavelength must be row or column vector')
        irradWavelength = irradWavelength(:);
        spectralIrradiance = p.Results.spectralIrradiance;
        assert(length(irradWavelength)==size(spectralIrradiance,1) ||...
            length(irradWavelength)==size(spectralIrradiance,2),...
            'irradWavelength and spectralIrradiance must be same length')
        if isrow(spectralIrradiance)
            spectralIrradiance = spectralIrradiance(:);
        elseif length(irradWavelength)==size(spectralIrradiance,2)
            spectralIrradiance = spectralIrradiance';
        end
    end
end

% check consistency of band designations
bandPass = getBands(p);

% interpolation function for spectral reflectance
DirectAndDiffuse = size(spectralReflectance,2)==2;
fitOptions = fitoptions('pchip');
if DirectAndDiffuse
    RdirF = fit(lambda,spectralReflectance(:,1),'pchip',fitOptions);
    RdifF = fit(lambda,spectralReflectance(:,2),'pchip',fitOptions);
else
    RdirF = fit(lambda,spectralReflectance,'pchip',fitOptions);
end

% interpolation function for spectral irradiance
if ~(ignoreSolar || (isequal(irradWavelength,Last_irradWavelength) &&...
        isequal(spectralIrradiance,Last_spectralIrradiance)))
    if size(spectralIrradiance,2)==2
        SolDirF = fit(irradWavelength,spectralIrradiance(:,1),'pchip',fitOptions);
        SolDifF = fit(irradWavelength,spectralIrradiance(:,2),'pchip',fitOptions);
    else
        SolDirF = fit(irradWavelength,spectralIrradiance,'pchip',fitOptions);
    end
    Last_irradWavelength = irradWavelength;
    Last_spectralIrradiance = spectralIrradiance;
end

% make sure bandPass values are within spectrum
assert(all(bandPass(:)>=min(lambda)),...
    'check input - some bandPass wavelengths < min(lambda) %g',min(lambda))
assert(all(bandPass(:)<=max(lambda)),...
    'check input - some bandPass wavelengths > max(lambda) %g',max(lambda))
assert(all(bandPass(:)>=min(irradWavelength)),...
    'check input - some bandPass wavelengths < min(irradWavelength) %g',...
    min(irradWavelength))
assert(all(bandPass(:)<=max(irradWavelength)),...
    'check input - some bandPass wavelengths > max(irradWavelength) %g',...
    max(irradWavelength))

R = zeros(size(bandPass,1),1);
for k=1:size(bandPass,1)
    if ignoreSolar
        R(k) = integral(@avgR,bandPass(k,1),bandPass(k,2))/...
            (bandPass(k,2)-bandPass(k,1));
    else
        denom = integral(@integrandS,bandPass(k,1),bandPass(k,2));
        if denom<=0
            R(k) = integral(@avgR,bandPass(k,1),bandPass(k,2))/...
                (bandPass(k,2)-bandPass(k,1));
        else
            num = integral(@integrandR,bandPass(k,1),bandPass(k,2));
            R(k) = num/denom;
        end
    end
end

if nargout>1
    varargout{1} = table(bandPass(:,1),bandPass(:,2),'VariableNames',...
        {'LowerWavelength','UpperWavelength'});
end

    function X = integrandR(wavelength)
        if DirectAndDiffuse
            X = RdirF(wavelength).*SolDirF(wavelength)+...
                RdifF(wavelength).*SolDifF(wavelength);
        else
            X = RdirF(wavelength).*SolDirF(wavelength);
        end
        X = X';
        X(X<0) = 0;
    end
    function S = integrandS(wavelength)
        if DirectAndDiffuse
            S = SolDirF(wavelength)+SolDifF(wavelength);
        else
            S = SolDirF(wavelength);
        end
        S = S';
        S(S<0) = 0;
    end
    function X = avgR(wavelength)
        X = RdirF(wavelength)';
    end
end

function bandPass = getBands(p)
% check consistence of band designations
assert(xor(isempty(p.Results.bandpass),isempty(p.Results.sensor)),...
    'either ''bandPass'' or ''sensor''/''bands'' must be specified, but not both')
if ~isempty(p.Results.bandpass)
    bandPass = p.Results.bandpass;
else
    T = SensorTable(p.Results.sensor,p.Results.lambdaUnit);
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