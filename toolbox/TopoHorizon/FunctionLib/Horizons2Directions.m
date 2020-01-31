function [sinF,sinB,azmF,azmB] = Horizons2Directions(planet,angToCalculate,lat,lon,Z)
% [sinF,sinB,azmF,azmB] = Horizons2Directions(planet,angToCalculate,lat,lon,Z)
%
% horizons in direction angToCalculate (e.g. angle toward sun)
% values returned are the sines of the horizon angles and the azimuths for
% each point (which can vary when dealing with projected or geolocated
% data)
%
% Input
%   planet - e.g. 'earth'
%   angToCalculate - nominal azimuth in forward direction, wrt grid (which
%       may not be wrt north)
%   lat, lon - latitude & longitude, degrees
%   Z - elevation, either single or double
%
% Output
%   sinF - grid of horizon sines, of size(Z), in forward direction
%   sinB - grid of horizon sines, of size(Z), in backward direction
%   azmF, azmB - forward and backward azimuths, of size(Z)
%       ("forward" is defined as the angle above which the points are
%       hidden when facing the sun at its azimuth)

[ hForward, hBackward, azmF, azmB ] =...
    horz2d( lat,lon,Z,angToCalculate,planet );
% sines of the angles
sinF = sind(hForward);
sinB = sind(hBackward);

end