function [ lat,lon] = MODtile2latlon( tile, varargin )
% [ lat,lon] = MODtile2latlon( tile, ['polygon' or 'makegrid'] )
%corner coordinates of MODIS tile, counter-clockwise from upper left
%
% Input
%   tile - MODIS tile in form 'hNNvNN'
% Optional input
%   'polygon' - if specified, close the lat-lon polygon
%   'makegrid' - instead return grids of latitude and longitude for each
%       pixel, must specify additional argument '250m' or '500m' or '1km'
%
% Output
%   lat, lon - if not makegrid option: vectors of latitudes and longitudes
%       of corners, in degrees, counter-clockwise from upper left
%       (Vectors of size(4,1), unless 'polygon' is specified, in which case
%       vectors are of size(5,1)
%   if 'makegrid' option is specified, then full grids of latitude and
%   longitude are returned, of size of the MODIS tile
%
% Note: truly the boundaries are slightly curved

% parse input
p = inputParser;
defaultResolution = '';
defaultPolygon = false;
addRequired(p,'tile',@ischar);
addParameter(p,'makegrid',defaultResolution,@ischar);
addParameter(p,'polygon',defaultPolygon,@islogical);
parse(p,tile,varargin{:});
resolution = p.Results.makegrid;
polygon = p.Results.polygon;

persistent already mstruct
if isempty(already)
    already = true;
    % hardwired, authalic radius corresponding to WGS84
    MODISgeoid = [6.371007181e+06 0];
    mstruct = defaultm('sinusoid');
    mstruct.geoid = MODISgeoid;
    mstruct = defaultm(mstruct);
end

% make a grid
if ~isempty(resolution)
    R = sinusoidProjMODtile(tile);
    switch resolution
        case '250m'
            RefMatrix = R.RefMatrix_250m;
            siz = [4800 4800];
        case '500m'
            RefMatrix = R.RefMatrix_500m;
            siz = [2400 2400];
        case '1km'
            RefMatrix = R.RefMatrix_1km;
            siz = [1200 1200];
        case '1000m'
            RefMatrix = R.RefMatrix_1km;
            siz = [1200 1200];
        otherwise
            error('''makegrid'', ''%s'' invalid, options are ''250m'', ''500m'' or ''1km''',...
                resolution)
    end
    % grid points in projected coordinates
    [x,y] = pixcenters(RefMatrix,siz,'makegrid');
    % convert to lat,lon
    [lat,lon] = projinv(mstruct,x,y);
    
else % just the corners of tile
    [ULx,ULy,tileHeight,tileWidth] = MODtile2xy(tile);
    x = [ULx ULx ULx+tileWidth ULx+tileWidth]';
    y = [ULy ULy-tileHeight ULy-tileHeight ULy]';
    
    % corresponding latitudes and longitudes
    [lat,lon] = projinv(mstruct,x,y);
    
    if polygon
        lat = cat(1,lat,lat(1));
        lon = cat(1,lon,lon(1));
    end
end
end