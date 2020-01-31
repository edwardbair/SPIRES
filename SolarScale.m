function [ values,varargout ] = SolarScale(varargin )
% [ values [,wavelength] ] = SolarScale([wavelength,]Name,Value,...)
%Returns scaled solar radiation such that the integral from 0.28 to 5.0 um is 1.0.
%Therefore, to get spectral distribution of TOA or surface irradiance,
%multiply the result from SolarScale by the integrated value (i.e., "solar
%constant" if you want TOA spectrum)
% OptionalInput
%   wavelength - scalar, vector, or matrix of desired wavelengths (if
%   omitted, uses all wavelengths in original data)
% Other input, name-value pairs in either order
%   'location' - either 'TOA' or 'surface' or 'ratio' (default surface)
%   'units' - wavelength units, no default to prevent errors
%
% Output
%   if 'TOA' or 'surface'
%       returns scaled value of solar radiation, so that integral from 0.28 to 5 um
%       is 1.0 (so one can multiply by a total to get a spectrum)
%   if 'ratio'
%       returns atmospheric transmittance as a function of wavelength
%
% Reference: http://rredc.nrel.gov/solar/spectra/am1.5/ASTMG173/ASTMG173.html,
% which is based on the SMARTS v 2.9.2 model
% Gueymard, C. A. (2004), The sun's total and spectral irradiance for solar
% energy applications and solar radiation models, Solar Energy, 76, 423-453,
% doi: 10.1016/j.solener.2003.08.039.
%
% Values from 4-5 um filled in with MODTRAN values and the ATRAN
% transmission estimates:
% http://rredc.nrel.gov/solar/spectra/am0/modtran.html
% https://atran.sofia.usra.edu/cgi-bin/atran/atran.cgi

persistent already FTOA Fsurface Fratio units TableWavelength

p = inputParser;
% inputs and defaults
defaultLocation = 'surface';
defaultUnits = '';
nonnegativeValidation = @(x) isnumeric(x) && all(x(:)>=0);
addOptional(p,'wavelength',[],nonnegativeValidation)
addParameter(p,'location',defaultLocation,@ischar)
addParameter(p,'units',defaultUnits,@ischar)
parse(p,varargin{:});

% check units
assert (~isempty(p.Results.units), '''units'' must be specified')

% set options
doSurface = strcmpi(p.Results.location,defaultLocation);
doRatio = strcmpi(p.Results.location,'ratio');
doTOA = strcmpi(p.Results.location,'toa');
assert(doSurface || doRatio || doTOA,...
    '''location'' must be either ''surface'' or ''TOA'' or ''ratio''')

% first pass or different units, read table and calculate functions
if isempty(already) || isempty(units) || ~strcmpi(units,p.Results.units)
    m = matfile('ASTMG173.mat'); % load the solar radiation data
    SolarIrrad = m.SolarIrrad;
    units = p.Results.units;
    TOA = SolarIrrad.TOA;
    sRad = SolarIrrad.surface;
    t = sRad>=TOA;
    ratio = sRad./TOA;
    minD = min(1-ratio(~t));
    TOA(t) = sRad(t)*(1+minD);
    if strcmpi(units,'nm')
        wv = SolarIrrad.wavelength;
    else
        wv = convertUnits(SolarIrrad.wavelength,'nm',units);
    end
    ftoa = fit(wv,TOA,'pchipinterp');
    fs = fit(wv,sRad,'pchipinterp');
    Fratio = griddedInterpolant(wv,sRad./TOA,'pchip','none');
    % scale TOA and surface irradiance to integrate to 1.0
    value = integrate(ftoa,max(wv),min(wv));
    FTOA = griddedInterpolant(wv,TOA/value,'pchip','none');
    value = integrate(fs,max(wv),min(wv));
    Fsurface = griddedInterpolant(wv,sRad/value,'pchip','none');
    TableWavelength = wv;
end

if isempty(p.Results.wavelength)
    wavelength = TableWavelength;
else
    wavelength = p.Results.wavelength;
end

% return values from interpolant
if doSurface
    values = Fsurface(wavelength);
elseif doTOA
    values = FTOA(wavelength);
elseif doRatio
    values = Fratio(wavelength);
else
    error('should not reach this statement, check code')
end
already = true;
% set values outside range to zero
values(isnan(values)) = 0;
% return wavelengths if requested
if nargout>1
    varargout{1} = wavelength;
end
end