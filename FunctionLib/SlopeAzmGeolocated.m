function [ slope,aspect ] = SlopeAzmGeolocated( lat, lon, Z, geoid )
% slopes and aspects of elevation grid in geolocated coordinates
%
%INPUT (matrices must be same size)
%   lat - latitude matrix, in degrees
%   lon - longitude matrix, in degrees
%   Z - elevation matrix (m)
%   geoid - vector of length 2 [SemimajorAxis Eccentricity]
%
%OUTPUT
%   slope - from horizontal, degrees
%   aspect - 0 is south, +counterclockwise, degrees

% check for sizes
if ~isequal(size(lat),size(lon),size(Z))
    error('lat lon & Z matrices must have same dimensions')
end

% slope & aspect
[aspect, slope, ~,~] = gradientm(lat,lon,Z,geoid);

% translate so that 0 is south, positive counter-clockwise
aspect = 180-aspect;

% flat areas have aspect 0, whereas gradientm returns NaN
aspect(isnan(aspect)) = 0;

end