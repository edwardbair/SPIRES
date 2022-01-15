function [ Structure ] = GridStructure( gridtype, gridsize, varargin)
% [ Structure ] = GridStructure( gridtype, gridsize, varargin)
%creates and fills template for variables for terrain structures
%(e.g. elevation, slope/aspect, horizons etc)
%
% INPUT
%   gridtype - either 'projected', 'geographic', or 'geolocated'
%   gridsize - size of grid (can be 3d but only first 2 used)
%
% OPTIONAL INPUT - name-value pairs (if omitted, empty ones are placed in
%   structure)
%   if 'projected'
%       'Projection', then projection structure
%       'RefMatrix', then 3x2 reference matrix
%   if 'geographic'
%       'RefMatrix', then 3x2 reference matrix
%       either . . .
%       'Geoid', then either 2-element vector [SemimajorAxis Eccentricity]
%       or 'Ellipsoid', then name of Ellipsoid recognized by referenceEllipsoid()
%   if 'geolocated'
%       either . . .
%       'Geoid', then either 2-element vector [SemimajorAxis Eccentricity]
%       or 'Ellipsoid', then name of Ellipsoid recognized by referenceEllipsoid()
%       'Latitude', matrix of latitudes
%       'Longitude', matrix of longitudes of same size
%
% OUTPUT
%   Structure - specific to input gridtype

optargs = nargin-2;
assert(optargs==0 || mod(optargs,2)==0,...
    'optional arguments must be name-value pairs')

% default values
proj = defaultm('pcarree');
proj = defaultm(proj);
refM = zeros(3,2);
latitude = [];
longitude = [];
E = referenceEllipsoid('wgs84');
geoid = [E.SemimajorAxis E.Eccentricity];

switch gridtype
    case 'projected'
        for k=1:2:optargs
            if strcmpi(varargin{k},'projection')
                proj = varargin{k+1};
            elseif strcmpi(varargin{k},'refmatrix')
                refM = varargin{k+1};
            else
                error('%s - unrecognized argument',varargin{k})
            end
        end
        RR = refmatToMapRasterReference(refM,gridsize);
        Structure = struct('gridtype',gridtype,...
            'Projection',proj,'RefMatrix',refM,'rasterReference',RR);
    case 'geographic'
        for k=1:2:optargs
            if strcmpi(varargin{k},'refmatrix')
                refM = varargin{k+1};
            elseif strcmpi(varargin{k},'geoid')
                geoid = varargin{k+1};
            elseif strncmpi(varargin{k},'ellips',6)
                E = referenceEllipsoid(varargin{k+1});
                geoid = [E.SemimajorAxis E.Eccentricity];
            else
                error('%s - unrecognized argument',varargin{k})
            end
        end
        RR = refmatToGeoRasterReference(refM,gridsize);
        Structure = struct('gridtype',gridtype,...
            'RefMatrix',refM,'rasterReference',RR,'Geoid',geoid);
    case 'geolocated'
        for k=1:2:optargs
            if strcmpi(varargin{k},'geoid')
                geoid = varargin{k+1};
            elseif strncmpi(varargin{k},'ellips',6)
                E = referenceEllipsoid(varargin{k+1});
                geoid = [E.SemimajorAxis E.Eccentricity];
            elseif strncmpi(varargin{k},'lat',3)
                latitude = varargin{k+1};
            elseif strncmpi(varargin{k},'lon',3)
                longitude = varargin{k+1};
            else
                error('%s - unrecognized argument',varargin{k})
            end
        end
        Structure = struct('gridtype',gridtype,...
            'Geoid',geoid,'Lat',latitude,'Lon',longitude);
    otherwise
        error('first argument must be ''projected'', ''geographic'', or ''geolocated''')
end

end