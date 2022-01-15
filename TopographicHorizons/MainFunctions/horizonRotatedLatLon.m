function [SForward,SBackward] = horizonRotatedLatLon(angToRotate,Z,R,useParallel,varargin)
% [SForward,SBackward] = horizonRotatedLatLon(angToRotate,Z,R,useParallel [,E])
%
% Horizons in direction after grid rotatated by angToRotate, for elevation
% grids in geographic format with uniform angular spacing in latitude and
% uniform angular spacing in longitude (but spacing for latitude and
% longitude can be different).
% Values returned are the horizon angles and the azimuths for each point
% (which can differ from some analytic solution based on rotation angle).
%
% To pick a rotation angle based on the desired azimuth, use the
% rotationAngleFromAzimuth function.
%
% Input
%   angToRotate - angle to rotate the grid
%   Z - elevation, either single or double
%   R - raster reference
%   useParallel - either true or false or 1 or 0, parallelizes the analyses
%       along the columns of the rotated elevation grid
% Optional input, if ellipsoid not in R.GeographicCRS
%   E - reference ellipsoid (output from referenceEllipsoid)
%
% Output, forward and backward structures with
%   sin of horizon angles
%   mean azimuth of angles
%   distances to horizons
%   ("forward" is defined as the angle above which the points are
%   hidden when facing the sun at its azimuth)

p = inputParser;
addRequired(p,'angToRotate',@(x) isnumeric(x) && isscalar(x) &&...
    abs(x)<=180)
addRequired(p,'Z',@(x) isnumeric(x) && ismatrix(x) && isfloat(x))
addRequired(p,'R',@(x) contains(class(x),'rasterref') &&...
    contains(class(x),'Geographic'))
addRequired(p,'useParallel',@(x) isscalar(x) &&...
    (islogical(x) || isnumeric(x)))
addOptional(p,'E',referenceEllipsoid('wgs84'),...
    @(x) contains(class(x),'referenceEllipsoid') ||...
    (isnumeric(x) && length(x)==2))

parse(p,angToRotate,Z,R,useParallel,varargin{:})
if isempty(p.Results.E) ||...
        (any(contains(fieldnames(R),'GeographicCRS')) &&...
        ~isempty(R.GeographicCRS))
    E = R.GeographicCRS.Spheroid;
else
    E = p.Results.E;
end

% latitudes and longitudes of the input grid
[xIntrinsic,yIntrinsic] = meshgrid(1:size(Z,2),1:size(Z,1));
[lat,lon] = intrinsicToGeographic(R,xIntrinsic,yIntrinsic);

if mod(angToRotate,90)==0
    [SForward,SBackward] = internal90(angToRotate,lat,lon,Z,E,...
        useParallel);
else
    [SForward,SBackward] = internalHorizonPts(angToRotate,lat,lon,Z,E,...
        useParallel);
end

end