function [mu0, phi0, airmass] = sunang( lat, lon, declin, omega, varargin )
% [mu0, phi0, airmass] = sunang( lat, lon, declin, omega, ... )
%Sun angles on horizontal surface, with optional correction for
%   refraction
%   (uses the distance function from the Mapping Toolbox)
%
%  Input (in degrees)
%   latitude
%   longitude
%   declin, solar declination
%   omega, longitude at which sun is vertical
%
%  Optional arguments that go after omega, in the following order:
%   zeroflag - if true, sets output negative cosines to zero and their
%       corresponding azimuths to NaN, default true
%       if false, returns negative cosines and their azimuths
%   (next optional arguments are a pair, supply neither or both)
%   P - pressure, to correct for atmospheric path length and
%       refraction, default ignore refraction
%   T - temperature in Kelvin, to correct for atmospheric path length and
%   	refraction, default ignore refraction
%  Output:
%   mu0, cosine of solar zenith angle
%   phi0, solar azimuth (degrees, from south, + ccw)
%   airmass - relative atmospheric path length, where 1.0 is the path
%       length at solar zenith angle = 0

% check for, and if necessary process, optional arguments
optargin=size(varargin,2);
assert( nargin-optargin >= 3,...
    'need 3 required arguments: latitude, declination, omega')
if optargin >= 1
    zeroflag = varargin{1};
    if optargin == 2
        error('if you specify pressure, you must specify temperature')
    elseif optargin == 3
        P = varargin{2};
        T = varargin{3};
    elseif optargin ~= 1
        error('too many optional arguments: can be 1 or 3')
    end
else
    % default is to set negative cosines to zero
    zeroflag = true;
end

if ~isscalar(lat)
    assert(isequal(size(lat),size(lon)),...
        'if not scalars, lat & lon must be same size')
    if ~isscalar(declin)
        assert(isequal(size(declin),size(omega),size(lat)),...
            'if not scalars, decline & omega must be same size as lat/lon')
    end
end

[arclen, phi0] = distance(lat, lon, declin, omega);
mu0 = cosd(arclen);
% translate so that 0 is south, positive counter-clockwise
phi0 = 180-phi0;

% atmospheric refraction
if optargin == 3
    mu0 = refracted(mu0, P/10, T);
end

% relative airmass
airmass = kasten(mu0);

% set negative cosines to zero
if zeroflag
    t = mu0 < 0;
    if nnz(t)
        mu0(t) = 0;
        phi0(t) = NaN;
    end
end

end