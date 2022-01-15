
function [ R,varargout ] = bandPassReflectance(reflTbl,solarTbl,varargin)
% [ R ] = bandPassReflectance(reflTbl,solarTbl,name/value pairs)
% [R, P] = bandPassReflectance(________________)
%integration over wavelength ranges from spectral reflectance and irradiance
%
%Inputs
% reflTbl - table with variable names 'wavelength' and 'reflectance', where
%   reflectance can have one or two columns: one column if reflectance is
%   total, or two columns with column 1 being direct reflectance at the
%   illumination angle acos(muS) and column 2 being diffuse reflectance
%   (units for wavelength must be specified in the
%   reflTbl.Properties.VariableUnits
% solarTbl - table or structure of solar radiation
%   empty, in which case the reflectance values are averaged across
%       wavelengths
%   table with 'wavelength' as one column and 'irradiance' as another either
%       as a single column or with two columns in which case the 2nd column
%       represents the diffuse irradiance;
%       in the table, the wavelengths can be different than those in the
%       reflTbl, but they must cover the range of the reflTbl and they must
%       have the same units as specified in
%       solarTbl.Properties.VariableUnits
%   structure with functions to return direct (on horizontal surface) and
%       diffuse irradiance with wavelength as the argument (in this case,
%       no checking for compatibility with units or range)
%   Note: the values of the solar radiation are used in a relative, not
%       absolute sense, so they can for example be normalized to the
%       maximum, normalized to integrate to 1.0 (as returned by SolarScale),
%       or they can represent geophysical values
%Optional input, name-value pairs (any non-ambiguous abbreviation works)
% either 'bandPass' or 'sensor'/'bands' must be specified but not both
% 'bandPass' - matrix of size Nx2, specifying the wavelength ranges of N
%   band passes, must be in the same units as lambdaUnit
% 'sensor' - followed by character string, any sensor from SensorTable.m
% 'bands' - followed by either numeric vector, cell vector, or categorical
%   vector of bands, or, if omitted, all bands for that sensor within the
%   range of lambda
%   (numeric vector can be used for those sensors whose bands are designated
%   just by number)
% 'mu0' - cosine of illumination angle on flat surface
% 'muS' - cosine of illumination angle on sloping surface
%Output
% R - band-averaged reflectance matrix, with N rows, the columns
%   representing each band
%Optional output
% P - prescription with functions for direct and diffuse irradiance, and
%   bandPass table with revised input bandPass values to include just the
%   bands that the wavelengths of the input reflectance spectra and
%   illumination spectra cover
%
%MATLAB requirements
% uses fit and integrate from the Curve Fitting Toolbox

% parse inputs
narginchk(5,13)
nargoutchk(0,2)
p = inputParser;
bpValidation = @(x) isnumeric(x) && all(x(:)>=0) && (size(x,2)==2 || size(x,1)==2);
bandValidation = @(x) isempty(x) || ((isrow(x) || iscolumn(x)) &&...
    ((isnumeric(x) && all(x(:)>0)) || iscell(x) || iscategorical(x)));
addRequired(p,'reflTbl',@istable)
addRequired(p,'solarTbl',@(x) istable(x) || isempty(x) || isstruct(x))
addParameter(p,validatestring('bandPass',{'bandp','bandPass'}),[],bpValidation)
addParameter(p,validatestring('sensor',{'sen','sensor'}),'',...
    @(x) ischar(x) || iscategorical(x));
addParameter(p,validatestring('bands',{'band','bands'}),[],bandValidation)
addParameter(p,'mus',0.5,@(x) isscalar(x) && x>=0 && x<=1)
addParameter(p,'mu0',0.5,@(x) isscalar(x) && x>=0 && x<=1)
parse(p,reflTbl,solarTbl,varargin{:})

muS = p.Results.mus;
mu0 = p.Results.mu0;
reflTbl = sortrows(p.Results.reflTbl);
if isempty(reflTbl.Properties.VariableUnits)
    warning('reflTbl doesn''t specify unit for wavelength, assuming ''um''')
    lambdaUnit = 'um';
else
    t = contains(reflTbl.Properties.VariableNames,'wavelength','ignoreCase',true);
    lambdaUnit = reflTbl.Properties.VariableUnits{t};
end

ignoreSolar = isempty(p.Results.solarTbl);
if ~ignoreSolar
    if isstruct(p.Results.solarTbl)
        Pstruct = p.Results.solarTbl;
        SolDirF = Pstruct.SolDirF;
        if isfield(Pstruct,'SolDifF')
            SolDifF = Pstruct.SolDifF;
            solarDiffuse = true;
        else
            solarDiffuse = false;
        end
        bandPass = Pstruct.bandPass;
        bandID = Pstruct.bandID;
    else % it's a table
        solarTbl = sortrows(p.Results.solarTbl);
        if isempty(solarTbl.Properties.VariableUnits)
            warning('solarTbl.Properties.VariableUnits doesn''t show wavelength units')
        else
            t = contains(solarTbl.Properties.VariableNames,'wavelength');
            solarWaveUnit = solarTbl.Properties.VariableUnits{t};
            %check units, allowing for 'mum' and 'um' to be the same
            assert(strcmpi(lambdaUnit,solarWaveUnit) ||...
                (contains(lambdaUnit,'um') && contains(solarWaveUnit,'um')),...
                'solar units %s not same as reflectance units %s',...
                solarWaveUnit,lambdaUnit)
        end
        % make sure solar table covers the wavelengths of reflTbl
        assert(min(solarTbl.wavelength)<=min(reflTbl.wavelength) &&...
            max(solarTbl.wavelength)>=max(reflTbl.wavelength),...
            'wavelengths of solar table (%g to %g) don''t cover wavelengths of reflectance table (%g tl %g)',...
            min(solarTbl.wavelength),max(solarTbl.wavelength),...
            min(reflTbl.wavelength),max(reflTbl.wavelength))
    end
end

% solar values
if ~ignoreSolar && ~exist('Pstruct','var')
    irradWavelength = double(solarTbl.wavelength);
    directIrradiance = double(solarTbl.irradiance(:,1));
    if size(solarTbl.irradiance,2)==2
        diffuseIrradiance = double(solarTbl.irradiance(:,2));
        solarDiffuse = any(diffuseIrradiance>0);
    else
        solarDiffuse = false;
    end
    SolDirF = fit(irradWavelength,directIrradiance,'pchip');
    if solarDiffuse
        SolDifF = fit(irradWavelength,diffuseIrradiance,'pchip');
    end
end

% reflectance values
reflWavelength = double(reflTbl.wavelength);
directReflectance = double(reflTbl.reflectance(:,1));
if size(reflTbl.reflectance,2)==2
    diffuseReflectance = double(reflTbl.reflectance(:,2));
    reflDiffuse = true;
else
    reflDiffuse = false;
end

%check if Diffuse, both diffuse irradiance and diffuse reflectance must be
%specified
assert(~xor(solarDiffuse,reflDiffuse),...
    'if diffuse reflectance is to be incorporated, both diffuse spectral reflectance and solar spectral diffuse irradiance must be provided')

%check consistency of band designations
if ~exist('Pstruct','var')
    [bandPass,bandID] = getBands(p);
end

% interpolation function(s) for irradiance and reflectance
if ~ignoreSolar
    upFluxDir = fit(reflWavelength,directReflectance.*SolDirF(reflWavelength),'pchip');
    if reflDiffuse
        upFluxDif = fit(reflWavelength,diffuseReflectance.*SolDifF(reflWavelength),'pchip');
    end
else
    ReflDirF = fit(reflWavelength,directRefl,'pchip');
end

bandReflectance = zeros(size(bandPass,1),1);
if ignoreSolar % i.e. integrate just the reflectance across the band passes
    bDiff = bandPass(:,2)-bandPass(:,1);
    for k=1:size(bandPass,1)
        bandReflectance(k) =...
            integrate(ReflDirF,bandPass(k,2),bandPass(k,1));
    end
    bandReflectance = bandReflectance./bDiff;
else
    TotalDirect = zeros(size(bandPass,1),1);
    TotalDiffuse = zeros(size(TotalDirect));
    for k=1:size(bandPass,1)
        TotalDirect(k) = integrate(SolDirF,bandPass(k,2),bandPass(k,1));
        if solarDiffuse
            TotalDiffuse(k) = integrate(SolDifF,bandPass(k,2),bandPass(k,1));
        end
    end
    
    if reflDiffuse
        for k=1:size(bandPass,1)
            bandReflectance(k) =...
                integrate(upFluxDir,bandPass(k,2),bandPass(k,1))*muS/mu0 +...
                integrate(upFluxDif,bandPass(k,2),bandPass(k,1));
        end
    else
        for k=1:size(bandPass,1)
            bandReflectance(k) =...
                integrate(upFluxDir,bandPass(k,2),bandPass(k,1))*muS/mu0;
        end
    end
    bandReflectance = bandReflectance./(TotalDirect+TotalDiffuse);
end

if nargout>1
    Pstruct.SolDirF = SolDirF;
    Pstruct.SolDifF = SolDifF;
    Pstruct.bandPass = bandPass;
    Pstruct.bandID = bandID;
    varargout{1} = Pstruct;
end
R = bandReflectance.';
end

function [bandPass,bandID] = getBands(p,lambdaUnit)
% check consistence of band designations
assert(xor(isempty(p.Results.bandPass),isempty(p.Results.sensor)),...
    'either ''bandPass'' or ''sensor''/''bands'' must be specified, but not both')
if ~isempty(p.Results.bandPass)
    bandPass = p.Results.bandPass;
    bandID = (1:size(bandPass,1))';
    % make sure bigger number is in column 2
    bandPass = sort(bandPass,2);
else
    T = SensorTable(p.Results.sensor,lambdaUnit);
    if isempty(p.Results.bands)
        % all bands
        bandPass = [T.LowerWavelength T.UpperWavelength];
        bandID = T.Band;
    else
        x = categorical(p.Results.bands);
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
t = bandPass(:,1)<min(p.Results.reflTbl.wavelength) |...
    bandPass(:,2)>max(p.Results.reflTbl.wavelength);
if any(t)
    warning('some bands outside range of reflectances, therefore truncated')
    bandPass = bandPass(~t,:);
    bandID = bandID(~t);
end
end