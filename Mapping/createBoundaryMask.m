function [ mask ] = createBoundaryMask(varargin )
% [ mask ] = createBoundaryMask(...)
%
%create logical mask for specified boundary
%
% This function is used to restrict processing in an image to a specified
% area.
%
% Input, all name-value pairs and all case-insensitive
%   'area' followed by name like 'Sierra' or 'HinduKush', or alternatively
%       by a shape structure with Geometry = 'Polygon' and either X,Y or
%       Lat,Lon coordinates. If the coordinates are X,Y, they are assumed
%       to be in the same projection as the mask to be created.
%   'proj' followed by a MATLAB projection structure, this can be omitted
%       if the mask is in geographic coordinates or if the projection is
%       the same as the input area
%   (in the following, either both 'RefMatrix' and 'imagesize', or
%   'rasterreference' alone)
%   'RefMatrix' followed by a 3x2 referencing matrix for the output mask
%       area
%   'size' followed by the size of the output mask
%   'rasterref' map or geographic raster reference object
%
% Output
%   mask - logical matrix of the right size, true within the boundary,
%       false outside

p = inputParser;
% validation functions
vRfcn =@(x) validateattributes(x,{'numeric'},{'size',[3 2]}); % ref matrix
vAfcn =@(x) ~isempty(x) && (isstruct(x) || ischar(x)); % area
vPfcn =@(x) isstruct(x); % projection structure
vZfcn =@(x) validateattributes(x,{'numeric'},{'size',[1 2]}); % image size
vRRfcn =@(x) isa(x,'map.rasterref.MapCellsReference') ||...
    isa(x,'map.rasterref.GeographicCellsReference'); % raster reference
% convert all to lower case
assert(~isempty(varargin),'no arguments, got to be some')
for k=1:length(varargin)
    if ischar(varargin{k})
        varargin{k} = lower(varargin{k});
    end
end
% defaults all empty (except 'area' which will be caught)
addParameter(p,'area',[],vAfcn);
addParameter(p,'proj',[],vPfcn);
addParameter(p,'refmatrix',[],vRfcn);
addParameter(p,'size',[],vZfcn);
addParameter(p,'rasterref',[],vRRfcn);

% parse inputs
parse(p,varargin{:})
assert(~isempty(p.Results.area),'''area'' must be specified')

if ischar(p.Results.area)
    % named area - truncate to 5 characters
    area = p.Results.area(1:5);
    m = matfile('BoundaryShapes.mat');
    switch area
        case 'sierr'
            B = m.SierraBoundary;
        case 'hindu'
            B = m.HinduKushBoundary;
        otherwise
            error('area ''%s'' not recognized',p.Results.area)
    end
else
    % shape structure, make sure valid
    B = p.Results.Area;
    assert(strcmp(B.Geometry,'Polygon'),...
        '''area'' if a structure must be a ''Polygon''')
    assert((isfield(B,'X') && isfield(B,'Y')) ||...
        (isfield(B,'Lat') && isfield(B,'Lon')),...
        '''area'' if a structure must have X-Y or Lat-Lon coordinates')
end
% projection, if specified, 'area' must be Lat-Lon
proj = p.Results.proj;
if ~isempty(proj)
    if isfield(B,'X')
        B.Lon = B.X;
        B.Lat = B.Y;
        B = rmfield(B,'X');
        B = rmfield(B,'Y');
        warning('''area'' structure has X-Y coordinates, assumed Lat-Lon')
    end
    % convert the boundary to the projection
    [xBoundary,yBoundary] = mfwdtran(proj,B.Lat,B.Lon);
end
% referencing matrix or raster reference
if ~isempty(p.Results.rasterref)
    rr = p.Results.rasterref;
    [xIntrinsic, yIntrinsic] = worldToIntrinsic(rr,xBoundary,yBoundary);
    mask = poly2mask(xIntrinsic,yIntrinsic,rr.RasterSize(1),rr.RasterSize(2));
elseif ~isempty(p.Results.refmatrix)
    RefMatrix = p.Results.refmatrix;
    rasterSize = p.Results.size;
    assert(~isempty(rasterSize),'if ''RefMatrix'' specified, so must ''size''')
    % pixel coordinates of the boundary
    [yIntrinsic,xIntrinsic] = map2pix(RefMatrix,xBoundary,yBoundary);
    mask = poly2mask(xIntrinsic,yIntrinsic,rasterSize(1),rasterSize(2));
else
    error(['either ''RefMatrix'' and ''size'' must be specified,'...
        ' or ''rasterref'' must be specified'])
end
% issue warning if both approaches specified
if ~isempty(p.Results.refmatrix) && ~isempty(p.Results.rasterref)
    warning(['both ''RefMatrix'' and ''rasterref'' specified, this is'...
        ' redundant, so ''rasterref'' used and ''RefMatrix'' ignored'])
end
end