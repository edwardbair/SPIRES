function [ R,varargout ] = bandPassRefl(lambda,lambdaUnit,spectralReflectance,solarRad,varargin)
% [ R ] = bandPassRefl(lambda, lambdaUnit, spectralReflectance,
%   spectralIrradiance, name/value pairs)
% [ R ] = bandPassRefl(lambda, lambdaUnit, spectralReflectance,
%   spectralIrradiance, name/value pairs)
% [R, bandPass] = bandPassRefl(________________)
%integration over wavelength ranges from spectral reflectance and irradiance
%
%Inputs
% lambda - wavelengths of spectralReflectance values, can be row or column
%   vector
% lambdaUnit - any length unit, typically 'um', 'mum', or 'nm'
% spectralReflectance - can be a vector of same size as lambda, or a matrix
%   of size [N length(lambda)] where the rows represent different
%   observations and the columns correspond to lambda
%   (to get separate values for direct and diffuse reflectance, specify
%   'diffuseReflectance' as an optional input and make spectralReflectance
%   correspond to direct)
% solarRad - can be one of the following:
%   empty, in which case the reflectance values are averaged across
%       wavelengths;
%   row or column vector of same length (L) as lambda or a matrix of size
%       [2 L] or [L 2], in which case the first row or column is the direct
%       irradiance and the second is diffuse;
%   table with 'wavelength' as one column and 'irradiance' as another either
%       as a single column or with two columns in which case the 2nd column
%       represents the diffuse irradiance;
%       in the table, the wavelengths can be different than lambda, but they
%       should cover the range [min(lambda) max(lambda)] and they must have
%       the same units as lambda
%   Note: the values of the solar radiation are used in a relative, not
%       absolute sense, so they can for example be normalized to the
%       maximum, or they can represent correct geophysical values
%Optional input, name-value pairs (any non-ambiguous abbreviation works)
% either 'bandPass' or 'sensor'/'bands' must be specified, but not both
% 'bandPass' - matrix of size Nx2, specifying the wavelength ranges of N
%   band passes, must be in the same units as lambdaUnit
% 'sensor' - followed by character string, any sensor from SensorTable.m
% 'band' - followed by either numeric vector, cell vector, or categorical
%   vector of bands, or, if omitted, all bands for that sensor within the
%   range of lambda
%   (numeric vector can be used for those sensors whose bands are designated just by number)
% 'diffuseReflectance' - for the same wavelengths as lambda, must be same
%   size(spectralReflectance)
% 'useParallel' - default false, but useful if many sets of spectral
%   reflectance values are to be considered
%
%Output
% R - band-averaged reflectance matrix, with N rows, the columns
%   representing each band
%Optional output
% bandPass - table with revised input bandPass values to include just the
%   bands that the wavelengths of the input reflectance spectra and
%   illumination spectra cover

% parse inputs
narginchk(6,14)
nargoutchk(0,2)
p = inputParser;
nonNegativeFcn = @(x) isnumeric(x)&&all(x(:)>=0);
bpValidation = @(x) isnumeric(x) && all(x(:)>=0) && (size(x,2)==2 || size(x,1)==2);
bandValidation = @(x) (isrow(x) || iscolumn(x)) &&...
    ((isnumeric(x) && all(x(:)>0)) || iscell(x) || iscategorical(x));
addRequired(p,'lambda',nonNegativeFcn)
addRequired(p,'lambdaUnit',@ischar)
addRequired(p,'spectralReflectance',nonNegativeFcn)
addRequired(p,'solarRad',@(x) (isnumeric(x) && all(x(:)>=0)) || istable(x))
addParameter(p,validatestring('bandpass',{'bandp','bandpass'}),[],bpValidation)
addParameter(p,validatestring('sensor',{'sen','sensor'}),'',...
    @(x) ischar(x) || iscategorical(x));
addParameter(p,'band',[],bandValidation)
addParameter(p,validatestring('diffusereflectance',...
    {'diffuser','dirffuserefl','diffusereflectance'}),[],nonNegativeFcn)
addParameter(p,'useparallel',false,@(x) isnumeric(x) || islogical(x))
parse(p,lambda,lambdaUnit,spectralReflectance,solarRad,varargin{:})

useParallel = logical(p.Results.useparallel);

% check sizes & make column vectors and Lx2 matrices
assert(isrow(p.Results.lambda) || iscolumn(p.Results.lambda),...
    'lambda must be row or column vector')
lambda = p.Results.lambda(:); %make it a column vector
assert(isequal(unique(lambda),sort(lambda)),...
    'input wavelengths must be in sorted ascending order')
spectralReflectance = p.Results.spectralReflectance;
if isvector(spectralReflectance)
    assert(length(lambda)==length(spectralReflectance),...
        'if spectralReflectance is a vector, it must the same length as lambda')
    if isrow(spectralReflectance)
        spectralReflectance = spectralReflectance(:); %column vector
    end
else % spectralReflectance is a matrix
    assert(size(spectralReflectance,2)==length(lambda),...
        'the number of columns in spectralReflectance must = length(lambda)')
    spectralReflectance = spectralReflectance'; %spectra down the columns
end
diffuseReflectance = p.Results.diffusereflectance;
if ~isempty(diffuseReflectance)
    assert(isequal(size(diffuseReflectance),size(spectralReflectance)),...
        'if specified, diffuse reflectance must have same size as spectral reflectance')
end

Diffuse = false;
Direct = false;
ignoreSolar = isempty(p.Results.solarRad);
if ignoreSolar
    commonWavelength = lambda(:);
    Direct = true;
elseif istable(p.Results.solarRad)
    solarTbl = sortrows(p.Results.solarRad);
    irradWavelength = solarTbl.wavelength;
    directIrradiance = solarTbl.irradiance(:,1);
    Direct = any(directIrradiance>0);
    if ~isvector(p.Results.solarRad.irradiance)
        diffuseIrradiance = solarTbl.irradiance(:,2);
        Diffuse = any(diffuseIrradiance>0);
    end
    commonWavelength = unique([lambda; irradWavelength]);
elseif ~isvector(p.Results.solarRad) %both direct and diffuse irradiance
    if size(p.Results.solarRad,1)==2 % irradiance along the rows
        directIrradiance = p.Results.solarRad(1,:)';
        diffuseIrradiance = p.Results.solarRad(2,:)';
    else
        directIrradiance = p.Results.solarRad(:,1);
        diffuseIrradiance = p.Results.solarRad(:,2);
    end
    assert(length(directIrradiance)==length(lambda),...
        'if specified as a matrix (instead of a table) length(solarRad) must = length(lambda)')
    Direct = any(directIrradiance>0);
    Diffuse = any(diffuseIrradiance>0);
    commonWavelength = lambda;
    irradWavelength = lambda;
elseif isvector(p.Results.solarRad)
    directIrradiance = p.Results.solarRad'; % column vector
    assert(length(directIrradiance)==length(lambda),...
        'if specified as a vector (instead of a table) length(solarRad) must = length(lambda)')
    Direct = any(directIrradiance>0);
    commonWavelength = lambda;
    irradWavelength = lambda;
else %shouldn't reach this statement
    assert(Direct || Diffuse,...
        'error: neither direct nor diffuse band pass reflectance is to be calculated')
end

%check if Diffuse, both diffuse irradiance and diffuse reflectance must be
%specified
if Diffuse && isempty(diffuseReflectance)
    Diffuse = false;
    warning ('if diffuse irradiance is provided, so must diffuse reflectance, so this diffuse irradiance is ignored')
end

% check consistency of band designations
[bandPass,bandID] = getBands(p);

% interpolation function(s) for irradiance
if ~ignoreSolar
    if Direct
        SolDirF = fit(irradWavelength,directIrradiance,'pchip');
    end
    if Diffuse
        SolDifF = fit(irradWavelength,diffuseIrradiance,'pchip');
    end
end

% go through all observations
if Direct
    nObs = size(spectralReflectance,2);
elseif Diffuse
    nObs = size(diffuseReflectance,2);
end
% initialize to NaN to bypass any NaN input
bandReflectance = nan(length(bandID),nObs);

if useParallel
    % interpolation function for spectral reflectance using parallel
    if ignoreSolar
        passVals = {lambda,bandPass};
        W = parallel.pool.Constant(passVals);
        parfor n=1:nObs
            bandReflectance(:,n) =...
                bandReflParallel(double(spectralReflectance(:,n)),W)
        end
    elseif Direct && Diffuse
        passVals = {lambda,commonWavelength,bandPass,SolDirF,SolDifF};
        W = parallel.pool.Constant(passVals);
        parfor n=1:nObs
            bandReflectance(:,n) = BothParallel(spectralReflectance(:,n),...
                diffuseReflectance(:,n),W);
        end
    elseif xor(Direct,Diffuse)
        if Direct
            reflMatrix = spectralReflectance;
            solF = SolDirF;
        else
            reflMatrix = diffuseReflectance;
            solF = SolDifF;
        end
        passVals = {lambda,commonWavelength,bandPass,solF};
        W = parallel.pool.Constant(passVals);
        parfor n=1:nObs
            bandReflectance(:,n) = JustOneParallel(reflMatrix(:,n),W);
        end
    else
        error('either Direct or Diffuse should have been set by here')
    end
else
    if ignoreSolar % i.e. integrate just the reflectance across the band passes
%         spectralReflectance = spectralReflectance'; % columns along the spectrum to be next in memory
        for n=1:nObs
            if ~any(isnan(spectralReflectance(:,n)))
                RspecF = fit(lambda(:),double(spectralReflectance(:,n)),'pchip');
                bDiff = (bandPass(:,2)-bandPass(:,1));
                for k=1:size(bandPass,1)
                    bandReflectance(k,n) =...
                        integrate(RspecF,bandPass(k,2),bandPass(k,1));
                end
                bandReflectance(:,n) = bandReflectance(:,n)./bDiff;
            end
        end
    elseif Direct && Diffuse
        for n=1:nObs
            bandReflectance(:,n) =...
                BothNonParallel(lambda,commonWavelength,bandPass,...
                double(spectralReflectance(:,n)),...
                double(diffuseReflectance(:,n)),SolDirF,SolDifF);
        end
    elseif xor(Direct,Diffuse)
        if Direct
            reflMatrix = spectralReflectance;
            solF = SolDirF;
        else
            reflMatrix = diffuseReflectance;
            solF = SolDifF;
        end
        for n=1:nObs
            bandReflectance(:,n) =...
                JustOneNonParallel(lambda,commonWavelength,bandPass,...
                reflMatrix(:,n),solF);
        end
    else
        error('either Direct or Diffuse must have been identified')
    end
end
if nargout>1
    varargout{1} = table(bandID, bandPass(:,1),bandPass(:,2),mean(bandPass,2),...
        'VariableNames',...
        {'Band','LowerWavelength','UpperWavelength','CentralWavelength'});
end
R = bandReflectance';
end

function [bandPass,bandID] = getBands(p)
% check consistence of band designations
assert(xor(isempty(p.Results.bandpass),isempty(p.Results.sensor)),...
    'either ''bandPass'' or ''sensor''/''bands'' must be specified, but not both')
if ~isempty(p.Results.bandpass)
    bandPass = p.Results.bandpass;
    bandID = (1:size(bandPass,1))';
    % make sure bigger number is in column 2
    bandPass = sort(bandPass,2);
else
    T = SensorTable(p.Results.sensor,p.Results.lambdaUnit);
    if isempty(p.Results.band)
        % all bands
        bandPass = [T.LowerWavelength T.UpperWavelength];
        bandID = T.Band;
    else
        x = categorical(p.Results.band);
        bandID = x(:)';
        bandPass = zeros(length(x),2);
        for k=1:length(bandID)
            b = find(T.Band==bandID(k));
            if isempty(b)
                warning(['band ' bandID(k) ' not found'])
            end
            bandPass(k,:) = [T.LowerWavelength(b) T.UpperWavelength(b)];
        end
    end
end
%eliminate the bands outside the wavelengths of the input
t = bandPass(:,1)<min(p.Results.lambda) | bandPass(:,2)>max(p.Results.lambda);
bandPass = bandPass(~t,:);
bandID = bandID(~t);
%also check if outside the range of the irradiance wavelengths
if ~isempty(p.Results.solarRad) && istable(p.Results.solarRad)
    t = bandPass(:,1)<min(p.Results.solarRad.wavelength) |...
        bandPass(:,2)>max(p.Results.solarRad.wavelength);
    bandPass = bandPass(~t,:);
    bandID = bandID(~t);
end
end

function R = BothNonParallel(wavelength,commonWavelength,bandPass,...
    directVector,diffuseVector,SolDirF,SolDifF)
R = nan(size(bandPass,1),1);
if ~any(isnan(directVector) | isnan(diffuseVector))
    RspecF = fit(wavelength,directVector,'pchip');
    RdifF = fit(wavelength,diffuseVector,'pchip');
    Rintg = fit(commonWavelength,...
        RspecF(commonWavelength).*SolDirF(commonWavelength)+...
        RdifF(commonWavelength).*SolDifF(commonWavelength),'pchip');
    for k=1:length(R)
        R(k) = integrate(Rintg,bandPass(k,2),bandPass(k,1))/...
            (integrate(SolDirF,bandPass(k,2),bandPass(k,1))+...
            integrate(SolDifF,bandPass(k,2),bandPass(k,1)));
    end
end
end

function R = JustOneNonParallel(wavelength,commonWavelength,bandPass,...
    reflVector,solF)
R = nan(size(bandPass,1),1);
if ~any(isnan(reflVector))
    Rrefl = fit(wavelength,reflVector,'pchip');
    Rintg = fit(commonWavelength,...
        Rrefl(commonWavelength).*solF(commonWavelength),'pchip');
    for k=1:length(R)
        R(k) = integrate(Rintg,bandPass(k,2),bandPass(k,1))/...
            integrate(solF,bandPass(k,2),bandPass(k,1));
    end
end
end

% function for the parallel option
function bandR = bandReflParallel(spectrum,W)
% bandR = bandRefl(spectrum,W)
% spectrum is one column of the reflectance matrix
% W contains the parallel pool constant
X = W.Value;
wavelength = X{1};
bandPass = X{2};
yfit = nan(size(bandPass,1),1);

%interpolating function
if ~any(isnan(spectrum))
    fR = fit(wavelength,spectrum,'pchip');
    bDiff = bandPass(:,2)-bandPass(:,1);
    for k=1:size(bandPass,1)
        yfit(k) = integrate(fR,bandPass(k,2),bandPass(k,1));
    end
    yfit = yfit./bDiff;
end
bandR = yfit;
end

function bandR = BothParallel(directRefl,diffuseRefl,W)
X = W.Value;
lambda = X{1};
commonW = X{2};
bandPass = X{3};
solDir = X{4};
solDif = X{5};
bandR = nan(size(bandPass,1),1);
if ~any(isnan(directRefl) | isnan(diffuseRefl))
    fDir = fit(lambda,directRefl,'pchip');
    fDif = fit(lambda,diffuseRefl,'pchip');
    Rintg = fit(commonW,fDir(commonW).*solDir(commonW)+...
        fDif(commonW).*solDif(commonW),'pchip');
    for k=1:size(bandPass,1)
        bandR(k) = integrate(Rintg,bandPass(k,2),bandPass(k,1))/...
            (integrate(solDir,bandPass(k,2),bandPass(k,1)))+...
            integrate(solDif,bandPass(k,2),bandPass(k,1));
    end
end
end

function bandR = JustOneParallel(reflVector,W)
X = W.Value;
lambda = X{1};
commonW = X{2};
bandPass = X{3};
solF = X{4};
bandR = nan(size(bandPass,1),1);
if ~any(isnan(reflVector))
    fRefl = fit(lambda,reflVector,'pchip');
    Rintg = fit(commonW,fRefl(commonW).*solF(commonW),'pchip');
    for k=1:size(bandPass,1)
        bandR(k) = integrate(Rintg,bandPass(k,2),bandPass(k,1))/...
            integrate(solF,bandPass(k,2),bandPass(k,1));
    end
end
end