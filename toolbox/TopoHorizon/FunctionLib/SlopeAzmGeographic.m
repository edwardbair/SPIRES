function [ slope, aspect ] = SlopeAzmGeographic( R, Z, varargin )
% slopes and aspects of elevation grid in geographic coordinates
%
%INPUT
%   R - referencing matrix for geographic grid, in degrees
%   Z - elevation grid (m)

%OPTIONAL INPUT specifies the ellipsoid
%   ellipsoid - either
%       a string identifying the ellipsoid (from MATLAB's
%           referenceEllipsoid function, with 'wgs84' the default)
%       or
%       a vector of length 2 [SemimajorAxis Eccentricity] - typically
%           this is in the 'geoid' field of the projection structure
%
%OUTPUT
%   slope - from horizontal, degrees
%   aspect - 0 is south, +counterclockwise, degrees

% lat, lon of each pixel
[lon, lat] = pixcenters(R,size(Z),'makegrid');

% geoid
% semimajor and semiminor axes
if isempty(varargin)
    E = referenceEllipsoid('wgs84');
    geoid = [E.SemimajorAxis E.Eccentricity];
else
    c = varargin;
    while iscell(c)
        c = c{1,1};
    end
    if isnumeric(c)
        geoid = c;
    elseif ischar(c)
        E = referenceEllipsoid(c);
        geoid = [E.SemimajorAxis E.Eccentricity];
    else
        celldisp(varargin)
        error('something wrong')
    end
end

[slope, aspect ] = SlopeAzmGeolocated(lat,lon,Z,geoid);

end