function [InRR,OutRR,planet,method,inLat,inLon,fillvalue ] =...
    parseInput(A,InR,InS,OutS,varargin)
%parse input values to produce output raster reference object

p = inputParser;
defaultPixelSize = NaN;
defaultXLimit = NaN;
defaultYLimit = defaultXLimit;
defaultAdjust = true;
defaultOrigin = 'ul';
defaultPlanet = 'earth';
defaultHemisphere = '';
defaultMethod = 'linear';
defaultRR = [];
defaultCells = true;
inLat = [];
inLon = [];

addRequired(p,'A',@(x) (isnumeric(x) || islogical(x) || iscategorical(x)) &&...
    (ismatrix(x) || ndims(x)==3))
addRequired(p,'InR',...
    @(x) contains(class(x),'rasterref') ||...
    (isnumeric(x) && isequal(size(x),[3 2])))
addRequired(p,'InS',@(x) isstruct(x) || isempty(x))
addRequired(p,'OutS',@(x) isstruct(x) || isempty(x))
addParameter(p,validatestring('pixelsize',{'pix','pixel','pixelsize'}),...
    defaultPixelSize,...
    @(x) isnumeric(x) && (isscalar(x) || length(x)==2))
addParameter(p,validatestring('XLimit',{'xli','xlim','xlimit'}),...
    defaultXLimit,@(x) isnumeric(x) && length(x)==2)
addParameter(p,validatestring('YLimit',{'yli','ylim','ylimit'}),...
    defaultYLimit,@(x) isnumeric(x) && length(x)==2)
addParameter(p,validatestring('adjust',{'adj','adjust'}),...
    defaultAdjust,@islogical)
addParameter(p,validatestring('cells',{'cel','cell','cells'}),...
    defaultCells,@islogical)
addParameter(p,validatestring('origin',{'ori','org','orig','origin'}),...
    defaultOrigin,@(x) ischar(x) &&...
    (strcmpi(x,'ul') || strcmpi(x,'ll') || strcmpi(x,'ur') || strcmpi(x,'lr')))
addParameter(p,validatestring('planet',{'pla','plan','planet'}),...
    defaultPlanet,@ischar)
addParameter(p,validatestring('hemisphere',{'hem','hemi','hemisphere'}),...
    defaultHemisphere,@ischar)
addParameter(p,validatestring('method',{'met','meth','method'}),...
    defaultMethod,@ischar)
addParameter(p,validatestring('rasterref',{'ras','rast','raster','rasterref'}),...
    defaultRR,@(x) contains(class(x),'rasterref'))
addParameter(p,validatestring('latitude',{'lat','lats','latitude'}),...
    [],@(x) isnumeric(x) && ismatrix(x))
addParameter(p,validatestring('longitude',{'lon','lons','longitude'}),...
    [],@(x) isnumeric(x) && ismatrix(x))
addParameter(p,validatestring('fillvalue',{'fil','fill','fillv','fillvalue'}),...
    [],@(x) isnumeric(x) && isscalar(x))


parse(p,A,InR,InS,OutS,varargin{:})

if isempty(p.Results.InS)
    InS = struct([]);
else
    InS = p.Results.InS;
end
if isempty(p.Results.OutS)
    OutS = struct([]);
else
    OutS = p.Results.OutS;
end

% which planet
planet = lower(p.Results.planet);

% interpolation method
if islogical(A)
    expectedMethods = {'nearest'};
elseif iscategorical(A)
    expectedMethods = {'nearest'};
elseif isempty(InR)
    expectedMethods = {'linear','nearest','natural'};
else
    expectedMethods = {'linear','nearest','cubic','spline','makima'};
end
method = validatestring(lower(p.Results.method),expectedMethods);

% bounding input coordinates
if strncmpi(p.Results.hemisphere,'nor',3)
    lonlim = [-180 180];
    latlim = [0 90];
elseif strncmpi(p.Results.hemisphere,'sou',3)
    lonlim = [-180 180];
    latlim = [-90 0];
elseif strncmpi(p.Results.hemisphere,'who',3) ||...
        strncmpi(p.Results.hemisphere,'bot',3)
    lonlim = [-180 180];
    latlim = [-90 90];
elseif isempty(p.Results.hemisphere)
    latlim = NaN;
    lonlim = NaN;
else
    error('''hemisphere'' ''%s'' not recognized',p.Results.hemisphere)
end

% check if size matches raster reference
sizeA = [size(p.Results.A,1) size(p.Results.A,2)];
if isempty(p.Results.InR)
    inLat = double(p.Results.latitude);
    inLon = double(p.Results.longitude);
    assert(isequal(size(inLat),sizeA) &&...
        isequal(size(inLon),sizeA),...
        'if input referencing matrix or raster reference is empty, latitude/longitude of same size as input raster must be specified')
    InRR = [];
    latlim = [min(inLat(:)) max(inLat(:))];
    lonlim = [min(inLon(:)) max(inLon(:))];
else
    assert(isempty(p.Results.latitude) && isempty(p.Results.longitude),...
        'if input referencing matrix or raster reference is specified, then input ''latitude'' and ''longitude'' must NOT be specified')
    if isempty(InS) % input is geographic, not projected
        if contains(class(p.Results.InR),'geographic','IgnoreCase',true)
            InRR = p.Results.InR;
        elseif ismatrix(p.Results.InR)
            if p.Results.cells
                InRR = refmatToGeoRasterReference(p.Results.InR,sizeA);
            else
                InRR = refmatToGeoRasterReference(p.Results.InR,sizeA,'postings');
            end
        elseif contains(class(p.Results.InR),'mapcells','IgnoreCase',true) ||...
                contains(class(p.Results.inR),'mappostings','IgnoreCase',true)
            error('if input A is geographic, input raster reference must be Geographic, not MapCells or MapPostings')
        end
        if isnan(latlim)
            [XIntrinsic,YIntrinsic] =...
                meshgrid([1 InRR.RasterSize(2)],[1 InRR.RasterSize(1)]);
            [latlim,lonlim] = intrinsicToGeographic(InRR,XIntrinsic,YIntrinsic);
            latlim = unique(latlim(:));
            lonlim = unique(lonlim(:));
        end
    else % input is projected
        if contains(class(p.Results.InR),'mapcells','IgnoreCase',true) ||...
                contains(class(p.Results.InR),'mappostings','IgnoreCase',true)
            InRR = p.Results.InR;
        elseif ismatrix(p.Results.InR)
            if p.Results.cells
                InRR = refmatToMapRasterReference(p.Results.InR,sizeA);
            else
                InRR = refmatToMapRasterReference(p.Results.InR,sizeA,'postings');
            end
        elseif contains(class(p.Results.InR),'geographic','IgnoreCase',true)
            error('if input A is in a projection, input raster reference must be Map, not Geographic')
        end
        if isnan(latlim)
            [XIntrinsic,YIntrinsic] =...
                meshgrid([1 InRR.RasterSize(2)],[1 InRR.RasterSize(1)]);
            [xWorld,yWorld] = intrinsicToWorld(InRR,XIntrinsic, YIntrinsic);
            try % minvtran fails on some projections if projection structure incorrect
                [latlim,lonlim] = minvtran(InS,xWorld,yWorld);
            catch % in those cases, use projinv instead, which also fails on some, like UTM
                [latlim,lonlim] = projinv(InS,xWorld,yWorld);
            end
        end
    end
end

% if output raster reference specified, other values default to it
% otherwise, parse other inputs
if isempty(p.Results.rasterref)
    % start of rows and cols
    if strncmpi(p.Results.origin,'l',1)
        startcol = 'south';
    else
        startcol = 'north';
    end
    if strcmpi(p.Results.origin(2),'l')
        startrow = 'west';
    else
        startrow = 'east';
    end
    
    % x- and y-limits
    if any(isnan(p.Results.xlimit(:))) || any(isnan(p.Results.ylimit(:)))
        if isempty(OutS)
            xlimit = [min(lonlim(:)) max(lonlim(:))];
            ylimit = [min(latlim(:)) max(latlim(:))];
        else
            try
                [x,y] = mfwdtran(OutS,latlim,lonlim);
            catch
                [x,y] = projfwd(OutS,latlim,lonlim);
            end
            xlimit = [min(x(:)) max(x(:))];
            ylimit = [min(y(:)) max(y(:))];
        end
    end
    if isnan(p.Results.xlimit)
        XLimit = xlimit;
    else
        XLimit = p.Results.xlimit;
    end
    if isnan(p.Results.ylimit)
        YLimit = ylimit;
    else
        YLimit = p.Results.ylimit;
    end
    
    % pixel size
    if isnan(p.Results.pixelsize)
        % pixel size is [height width], default is preserve image size
        pixelsize(1) = (YLimit(2)-YLimit(1))/size(A,1);
        pixelsize(2) = (XLimit(2)-XLimit(1))/size(A,2);
    else
        pixelsize = p.Results.pixelsize;
        if isscalar(pixelsize)
            pixelsize(2) = pixelsize(1);
        end
    end
    
    % adjust limits to multiple of pixel size
    if p.Results.adjust
        [XLimit,YLimit] = adjustLimits(XLimit,YLimit,pixelsize);
    end
    nRows = (YLimit(2)-YLimit(1))/pixelsize(1);
    nCols = (XLimit(2)-XLimit(1))/pixelsize(2);
    if nRows~=round(nRows)
        YLimit(2) = YLimit(1)+round(nRows)*pixelsize(1);
        nRows = round((YLimit(2)-YLimit(1))/pixelsize(1));
    end
    if nCols~=round(nCols)
        XLimit(2) = XLimit(1)+round(nCols)*pixelsize(2);
        nCols = round((XLimit(2)-XLimit(1))/pixelsize(2));
    end
    
    % output raster references
    if isempty(OutS)
        if p.Results.cells
            OutRR = georefcells(YLimit,XLimit,[nRows nCols],...
                'ColumnsStartFrom',startcol,'RowsStartFrom',startrow);
        else
            OutRR = georefpostings(YLimit,XLimit,[nRows nCols],...
                'ColumnsStartFrom',startcol,'RowsStartFrom',startrow);
        end
    else
        if p.Results.cells
            OutRR = maprefcells(XLimit,YLimit,[nRows nCols],...
                'ColumnsStartFrom',startcol,'RowsStartFrom',startrow);
        else
            OutRR = maprefpostings(XLimit,YLimit,[nRows nCols],...
                'ColumnsStartFrom',startcol,'RowsStartFrom',startrow);
        end
    end
else
    OutRR = p.Results.rasterref;
end

% fill value
fillvalue = p.Results.fillvalue;

end