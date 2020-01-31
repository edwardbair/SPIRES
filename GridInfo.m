function [ Structure ] = GridInfo( gridtype, varargin)
% [ Structure ] = GridInfo( gridtype, varargin)
%creates and fills template for variables for terrain structures
%(e.g. elevation, slope/aspect, horizons etc)
%
% INPUT
%   gridtype - either 'projected', 'geographic', or 'geolocated'
%
% OPTIONAL INPUT - name-value pairs, case-insensitive, generally any
%       unambiguous abbrevation of 3 or more letters works
%   'gridsize' then size of grid as vector of length 2 (will convert
%       RefMatrix to rasterReference if provided)
%   'Ellipsoid', as a referenceEllipsoid class, defaults to WGS84, other
%       ellipsoids or even other planets can be used (see documentation on
%       referenceEllipsoid)
%   if 'projected'
%       'Projection', then projection structure
%       'RefMatrix', then 3x2 reference matrix of projection coordinates
%       or 'rasterReference', in which case a variable of class
%           map.rasterref.MapCellsReference or map.rasterref.MapPostingsReference
%       (usually either 'RefMatrix' or 'rasterReference' is provided, not
%       both)
%   if 'geographic'
%       'RefMatrix', then 3x2 reference matrix of latitudes and longitudes
%       or 'rasterReference', in which case a variable of class
%           map.rasterref.GeographicCellsReference or
%           map.rasterref.GeographicPostingsReference
%       (usually either 'RefMatrix' or 'rasterReference' is provided, not
%       both)
%   if 'geolocated'
%       'Latitude', matrix of latitudes of same size as elevation grid
%       'Longitude', matrix of longitudes of same size as elevation grid
%
% OUTPUT
%   Structure - specific to input gridtype

optargs = nargin-1;
assert(optargs==0 || mod(optargs,2)==0,...
    'optional arguments must be name-value pairs')

%inputs
defaultE = referenceEllipsoid('wgs84');
p = inputParser;
addRequired(p,'gridtype',@ischar)
addParameter(p,validatestring('projection',{'pr','proj','projection'}),'',@isstruct)
addParameter(p,validatestring('refmatrix',{'ref','refm','refmatrix'}),...
    [],@(x) isnumeric(x) && isequal(size(x),[3 2]))
addParameter(p,validatestring('rasterreference',{'ras','rast','rasterref','rasterreference'}),...
    [],@(x) contains(class(x),'rasterref','IgnoreCase',true))
addParameter(p,validatestring('ellipsoid',{'ell','ellip','ellipsoid'}),...
    defaultE,@(x) contains(class(x),'referenceEllipsoid','IgnoreCase',true))
addParameter(p,validatestring('latitude',{'lat','latitude'}),...
    [],@ismatrix)
addParameter(p,validatestring('longitude',{'lon','long','longitude'}),...
    [],@ismatrix)
parse(p,gridtype,varargin{:})
gridtype = p.Results.gridtype;

% warning if both RefMatrix and rasterReference are specified
if ~isempty(p.Results.refmatrix) && ~isempty(p.Results.rasterreference)
    warning('''RefMatrix'' and ''rasterReference'' are both specified, generally you just specify one')
end

RR = [];
RefMatrix = [];
switch gridtype
    case 'projected'
        assert(~isempty(p.Results.refmatrix) || ~isempty(p.Results.rasterreference),...
            'if gridtype is ''projected'', ''RefMatrix'' or ''rasterReference'' must be specified')
        if ~isempty(p.Results.refmatrix)
            RefMatrix = p.Results.refmatrix;
            if ~isempty(p.Results.gridsize)
                RR = refmatToMapRasterReference(RefMatrix,p.Results.gridsize(:)');
            end
        end
        if ~isempty(p.Results.rasterreference)
            RR = p.Results.rasterreference;
            assert(contains(class(RR),'map','IgnoreCase',true),...'
                'if ''rasterReference'' specified for ''projected'' gridtype, must be class MapRasterReference, not %s',...
                class(RR))
            if isempty(RefMatrix)
                RefMatrix = RasterRef2RefMat(RR);
            end
        end
        Structure = struct('gridtype',gridtype,'Projection',p.Results.projection,...
            'RefMatrix',RefMatrix);
        if ~isempty(RR)
            Structure.rasterReference = RR;
        end
    case 'geographic'
        assert(~isempty(p.Results.refmatrix) || ~isempty(p.Results.rasterreference),...
            'if gridtype is ''geographic'', referencing matrix or rasterReference must be specified')
        if ~isempty(p.Results.refmatrix)
            RefMatrix = p.Results.refmatrix;
            if isempty(p.Results.gridsize)
                RR = [];
            else
                RR = refmatToGeoRasterReference(RefMatrix,p.Results.gridsize(:)');
            end
        end
        if ~isempty(p.Results.rasterreference)
            RR = p.Results.rasterreference;
            assert(contains(class(RR),'geo','IgnoreCase',true),...'
                'if ''rasterReference'' specified for ''geographic'' gridtype, must be class GeographicRasterReference, not %s',...
                class(RR))
            if isempty(RefMatrix)
                RefMatrix = RasterRef2RefMatrix(RR);
            end
        end
        Structure = struct('gridtype',gridtype,'Ellipsoid',p.Results.ellipsoid,...
            'RefMatrix',RefMatrix);
        if ~isempty(RR)
            Structure.rasterReference = RR;
        end
    case 'geolocated'
        Latitude = p.Results.latitude;
        Longitude = p.Results.longitude;
        assert(~isempty(Latitude) && ~isempty(Longitude),...
            'for ''geolocated'' gridtype, ''Latitude'' and ''Longtiude'' grids must be provided')
        assert(isequal(size(Latitude),size(Longitude)),...
            'for ''geolocated'' gridtype, ''Latitude'' and ''Longtiude'' grids must have same size')
        Structure = struct('gridtype',gridtype,...
            'Ellipsoid',p.Results.ellipsoid,'Latitude',Latitude,'Longitude',Longitude);
    otherwise
        error('first argument, gridtype, must be ''projected'', ''geographic'', or ''geolocated''')
end

end