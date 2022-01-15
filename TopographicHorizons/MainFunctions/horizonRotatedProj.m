function [SForward,SBackward] = horizonRotatedProj(angToRotate,Z,R,useParallel,varargin)
% [SForward,SBackward] = horizonRotatedProj(angToRotate,Z,R,useParallel [,proj])
%
% Horizons in direction after grid rotatated by angToRotate, for elevation
% grids in a projection.
% Values returned are the horizon angles and the azimuths for each point
% (which can differ from some analytic solution based on rotation angle
% depending on projection).
%
% To pick a rotation angle based on the desired azimuth, use the
% rotationAngleFromAzimuth function
%
% Input
%   angToRotate - angle to rotate the grid
%   Z - elevation, either single or double
%   R - raster reference
%   useParallel - either true or false or 1 or 0, parallelizes the analyses
%       along the columns of the rotated elevation grid, default false
% Optional input (if not specified in R.ProjectedCRS)
%   proj - projection structure or projcrs object
%
% Output
%   Structures for the forward and backward directions, with fields
%   including:
%   elevation angles in degrees of forward and backward horizons
%   distances to the horizon for each point
%   azimuths to horizons, averaged over the grid

p = inputParser;
addRequired(p,'angToRotate',@(x) isnumeric(x) && isscalar(x) &&...
    abs(x)<=180)
addRequired(p,'Z',@(x) isnumeric(x) && ismatrix(x) && isfloat(x))
addRequired(p,'R',@(x) contains(class(x),'rasterref') &&...
    contains(class(x),'Map'))
addRequired(p,'useParallel',@(x) isscalar(x) &&...
    (islogical(x) || isnumeric(x)))
addOptional(p,'proj',struct([]),@(x) contains(class(x),'projcrs') || isstruct(x))
parse(p,angToRotate,Z,R,useParallel,varargin{:})
if isempty(p.Results.proj)
    proj = R.ProjectedCRS;
    assert(~isempty(proj),'raster reference lacks a valid ProjectedCRS')
else
    proj = p.Results.proj;
end

% latitudes and longitudes of the input grid
[xIntrinsic,yIntrinsic] = meshgrid(1:R.RasterSize(2),1:R.RasterSize(1));
[xW,yW] = intrinsicToWorld(R,xIntrinsic,yIntrinsic);
[lat,lon] = projinv(proj,xW,yW);

if isstruct(proj)
    E = proj.geoid;
    % E = [SemimajorAxis Eccentricity]
    refEllipsoid = oblateSpheroid;
    refEllipsoid.SemimajorAxis = E(1);
    refEllipsoid.Eccentricity = E(2);
else
    if isempty(proj.GeographicCRS)
        refEllipsoid = referenceEllipsoid('earth');
    else
        refEllipsoid = proj.GeographicCRS.Spheroid;
    end
end

if mod(angToRotate,90)==0
    [SForward,SBackward] = internal90(angToRotate,lat,lon,Z,...
        refEllipsoid,useParallel);
else
    [SForward,SBackward] = internalHorizonPts(angToRotate,lat,lon,Z,...
        refEllipsoid,useParallel);
end

end