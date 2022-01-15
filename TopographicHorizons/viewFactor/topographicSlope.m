function [slope,aspect,varargout] = topographicSlope(Z,R,varargin)
% [slope,aspect] = topographicSlope(Z,R)
% [slope,aspect] = topographicSlope(Z,R,Name-Value pairs)
% [slope,aspect,gradN,gradE] = topographicSlope(Z,R,...)
%
%calculates slopes and aspects for elevation grid in all azimuths
%
%Input
%   Z - elevation grid in meters, single or double precision
%   R - raster reference for the elevation grid, can be a MapCellsReference
%       or MapPostingsReference object, in which case the grid is projected
%       and distances between points are based on the projected coordinates,
%       or a GeographicCellsReference or GeographicPostingsReference object,
%       in which case the distances between points are calculated from the
%       great-circle distances
%Optional input, name-value pairs, case-insensitive
%   'proj' - projection structure or projcrs object, needed only if R is a
%       MapCellsReference or MapPostingsReference object and its
%       ProjectedCRS field is missing or empty
%   'planet' - applicable only if R is a GeographicCellsReference or
%       GeographicPostingsReference object whose GeographicCRS field is
%       missing or empty, default 'wgs84', any ellipsoid
%       or planet recognized by referenceEllipsoid is recognized
%Output
%   slope - matrix of slopes corresponding to Z, in degrees upward from
%       horizontal
%   aspect - matrix of downhill directions, either ±180° from south,
%       positive east and negative west, or 0° to 360° clockwise from
%       north, depending on how the function azimuthPreference is set
%Optional output
%   gradN - north components of the gradient
%   gradE - east components of the gradient

% defaults
defaultPlanet = 'wgs84';
minargs = 2;
maxargs = 6;
narginchk(minargs,maxargs)

p = inputParser;
addRequired(p,'Z',@(x) isnumeric(x) && ismatrix(x))
addRequired(p,'R',@(x) contains(class(x),'rasterref'))
addParameter(p,'planet',defaultPlanet,@ischar)
addParameter(p,'proj',struct([]),@(x) isstruct(x) || contains(class(x),'projcrs'))
parse(p,Z,R,varargin{:});

planet = p.Results.planet;
if contains(class(R),'geographic','IgnoreCase',true)
    if ~isempty(p.Results.proj)
        warning('R (2nd argument) is of class %s, so ''proj'' should be empty, it''s ignored',class(R))
    end
    proj = struct([]);
else
    if any(contains(fieldnames(R),'ProjectedCRS')) && ~isempty(R.ProjectedCRS)
        proj = R.ProjectedCRS;
    else
        assert(~isempty(p.Results.proj),...
            'R (2nd argument) is of class %s and doesn''t contain a ProjectedCRS field, so ''proj'' must be specified',class(R))
        proj = p.Results.proj;
    end
end

useLatLon = contains(class(R),'geographic','IgnoreCase',true);
if useLatLon
    % check if R contains a GeographicCRS, which contains the ellipsoid
    if any(contains(fieldnames(R),'GeographicCRS')) && ~isempty(R.GeographicCRS)
        [aspect,slope,gradN,gradE] = gradientm(Z,R);
    else
        E = referenceEllipsoid(planet);
        % slope and aspect
        [aspect,slope,gradN,gradE] = gradientm(Z,R,E);
    end
else
    % latitudes and longitudes of the input grid
    [xIntrinsic,yIntrinsic] = meshgrid(1:R.RasterSize(2),1:R.RasterSize(1));
    [xW,yW] = intrinsicToWorld(R,xIntrinsic,yIntrinsic);
    % check if R contains a ProjectedCRS, which contains the projection and
    % the ellipsoid
    if any(contains(fieldnames(R),'ProjectedCRS')) && ~isempty(R.ProjectedCRS)
        refEllipsoid = R.ProjectedCRS.GeographicCRS.Spheroid;
    elseif contains(class(proj),'projcrs')
        refEllipsoid = proj.GeographicCRS.Spheroid;
    else % proj is a projection structure
        assert(isstruct(proj),'expecting ''proj'' to be a projection structure')
        E = proj.geoid;
        % E = [SemimajorAxis Eccentricity]
        refEllipsoid = oblateSpheroid;
        refEllipsoid.SemimajorAxis = E(1);
        refEllipsoid.Eccentricity = E(2);
    end
    % try alternatives to get lat/lon because of issues in versions of
    % MATLAB earlier than 2020b
    try
        [lat,lon] = projinv(proj,xW,yW);
    catch
        [lat,lon] = minvtran(proj,xW,yW); %#ok<MINVT>
    end
    [aspect,slope,gradN,gradE] = gradientm(lat,lon,Z,refEllipsoid);
end
% if azimuthPrefernce, translate so that 0 is south, positive counter-clockwise
if azimuthPreference
    aspect = 180-aspect;
end
% set flat areas aspect to zero, whereas gradientm returns NaN
aspect(isnan(aspect)) = 0;

% northern and eastern components of gradient
if nargout>2
    varargout{1} = gradN;
    if nargout>3
        varargout{2} = gradE;
    end
end
end