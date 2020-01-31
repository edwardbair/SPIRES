function [ slope,aspect ] = SlopeAzmProjected( R, P, Z )
% slopes and aspects of elevation grid in specified projection
%
%INPUT
%   R - referencing matrix for grid
%   P - projection structure
%   Z - elevation grid, oriented toward north
%
%OUTPUT
%   slope - from horizontal, degrees
%   aspect - 0 is south, +counterclockwise, 

% geographic coordinates
[x,y] = pixcenters(R,size(Z),'makegrid');
[lat,lon] = projinv(P,x,y);

[slope, aspect ] = SlopeAzmGeolocated(lat,lon,Z,P.geoid);

end