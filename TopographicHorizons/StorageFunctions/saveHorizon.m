function saveHorizon(filename,R,azm,horizons,varargin)
% saveHorizon(filename,R,projOrGeoid,azm,horizons,varargin)
%save horizons in format specified by filename:
%options are geotiff, HDF 5, or MATLAB
%if geotiff, only option is to read the whole file but then can be
%   processed by other software
%if HDF 5 (.h5), output is stored in block compressed form so horizons for
%   a specific azimuth can be retrieved without reading the whole file
%if MATLAB (.mat), output is stored as an interpolating function so that
%   horizons can be retrieved for specific azimuths and locations
%
%input
% filename, either .tif, .h5, or .mat, include full path if needed (if .tif
%   or .mat and distance output is specified, two output files will be
%   created, one filename containing 'horizon' and the other 'distance')
% R - raster reference for the horizon grid, can be a MapCellsReference
%   or MapPostingsReference object, in which case the grid the projection
%   argument is needed, or a GeographicCellsReference or
%   GeographicPostingsReference object, in which case the geoid can be
%   provided (default WGS84)
%   (If it's a 'postings' object, will be converted to 'cells')
% azm - vector of azimuths of horizons in degrees either from -180 to +180
%   ccw from north, with 0 south, or cw from 0 to 360 with 0 north, depending
%   on how the azimuthPreference function is set (output from horizonPostProcess
%   or raw values from horizonAllDirections)
% horizons - horizons in directions azm as BSQ array (output from
%   horizonAllDirections), one x-y or lat-lon plane for each azimuth, as
%   floating point numbers
%optional input, name/value pairs (unambiguous abbreviations of 4 or more
%   characters work)
% 'distances' - distances in directions azm as BSQ file (output from
%   horizonAllDirections), one x-y or lat-lon plane for each azimuth, as
%   floating point numbers, must be same size as horizons array
% 'proj' - projection structure or projCRS object, needed if R is a Map
%   Reference without a projcrs field
% 'geoid' - geoid as a name (e.g. 'wgs84'), or as a referenceEllipsoid
%   (output from referenceEllipsoid) or a 2-length vector of form
%   [SemimajorAxis Eccentricity]
% 'GeoKeyDirectoryTag' if output is '.tif', use the information from the
%   input .tif file (I=geotiffinfo(file), see I.GeoTIFFTags.GeoKeyDirectoryTag)
%

p = inputParser;
addRequired(p,'filename',@ischar)
addRequired(p,'R',@(x) contains(class(x),'rasterref'))
addRequired(p,'azm',@(x) isnumeric(x) && isvector(x))
addRequired(p,'horizons',@(x) isnumeric(x) && ndims(x)==3)
addParameter(p,validatestring('distances',{'dist','distance','distances'}),...
    [],@(x) isnumeric(x) && ndims(x)==3)
addParameter(p,validatestring('projection',{'proj','projection'}),...
    struct([]),@isstruct)
addParameter(p,validatestring('geoid',{'geoi','geoid'}),[],...
    @(x) isempty(x) || (isnumeric(x) && length(x)==2) || ischar(x) ||...
    contains(class(x),'referenceEllipsoid'))
addParameter(p,validatestring('geokeydirectorytag',{'geok','geokeydirectorytag'}),...
    struct([]),@isstruct);
parse(p,filename,R,azm,horizons,varargin{:})

% file must be .tif, .h5, or .mat
[~,~,ext] = fileparts(filename);
geotiff = false;
hdf5 = false;
mat = false;
switch ext
    case '.tif'
        geotiff = true;
        GeoKeyDirectoryTag = p.Results.geokeydirectorytag;
        GeoKeyDirectoryTag.GTRasterTypeGeoKey = 1;
    case '.h5'
        hdf5 = true;
    case '.mat'
        mat = true;
    otherwise
        error('filename must be .tif, .h5, or .mat')
end
% distances must be same size if entered
if ~isempty(p.Results.distances)
    distances = p.Results.distances;
    assert(isequal(size(distances),size(horizons)),...
        'if specified, distances must be same size as horizons')
end
% change R to 'cells' if currently 'postings
R = postings2cells(R);
% and identify whether projection or geoid is needed
if contains(class(R),'Map') && contains(class(R),'Reference')
    % map so must have projection
    if ~geotiff
        assert(isstruct(p.Results.projection),...
            'input rasterref is a % so ''projection'' argument is needed',class(R));
        proj = p.Results.projection;
    end
elseif contains(class(R),'geographic','IgnoreCase',true)
    % geographic so must have a geoid
    if ~geotiff
        if ischar(p.Results.geoid)
            E = referenceEllipsoid(p.Results.geoid);
            geoid = [E.SemimajorAxis E.Eccentricity];
        elseif isnumeric(p.Results.geoid)
            geoid = p.Results.geoid;
        elseif contains(class(p.Results.geoid),'referenceEllipsoid','IgnoreCase',true)
            E = p.Results.geoid;
            geoid = [E.SemimajorAxis E.Eccentricity];
        end
    end
else %shouldn't reach
    error('3rd argument unexpected type')
end

assert(isfloat(horizons),'input horizon array should be floating point (I could revise the code)')
if exist('distances','var')
    distances = single(distances);
end

if hdf5
    writeHDF5;
elseif geotiff
    writegeotiff;
elseif mat
    writeMAT;
end

    function writeMAT()
        % file name should include 'horizon'
        if contains(filename,'horizon','IgnoreCase',true) ||...
                contains(filename,'horz','IgnoreCase',true) ||...
                contains(filename,'hor','IgnoreCase',true)
            fn = filename;
        else % add 'horizon' to filename
            [fpath,fname,xt] = fileparts(filename);
            fname = [fname '_horizon'];
            fn = fullfile(fpath,[fname xt]);
        end
        %create lookup function
        [r,c,a] = ndgrid(1:size(horizons,1),1:size(horizons,2),azm);
        Fhorizon = griddedInterpolant(r,c,a,horizons,'linear','nearest');
        % write file
        writeProjection = false;
        writeGeoid = false;
        if exist('proj','var')
            projection = proj;
            writeProjection = true;
        end
        if exist('geoid','var')
            writeGeoid = true;
        end
        if writeProjection
            save(fn,'-v7.3','R','projection','Fhorizon')
        elseif writeGeoid
            save(fn,'-v7.3','R','geoid','Fhorizon')
        else
            error('either ''projection'' or ''geoid'' must be specified')
        end
        
        % same for distances
        if exist('distances','var')
            % change filename to include 'distance'
            [fpath,fname,xt] = fileparts(fn);
            if contains(fname,'horizon','IgnoreCase',true)
                fname = replace(lower(fname),'horizon','');
            elseif contains(lower(fname),'horz','IgnoreCase',true)
                fname = replace(lower(fname),'horz','');
            elseif contains(fname,'hor','IgnoreCase',true)
                fname = replace(lower(fname),'hor','');
            end
            fname = [fname '_distance'];
            fname = replace(fname,'__','_');
            fn = fullfile(fpath,[fname xt]);
            % interpolating function
            Fdistance = griddedInterpolant(r,c,a,single(distances),'linear','nearest');
            if writeProjection
                save(fn,'-v7.3','R','projection','Fdistance')
            elseif writeGeoid
                save(fn,'-v7.3','R','geoid','Fdistance')
            else
                error('either ''projection'' or ''geoid'' must be specified')
            end
        end
    end

    function writegeotiff()
        % file name should include 'horizon'
        if contains(filename,'horizon','IgnoreCase',true) ||...
                contains(filename,'horz','IgnoreCase',true) ||...
                contains(filename,'hor','IgnoreCase',true)
            fn = filename;
        else % add 'horizon' to filename
            [fpath,fname,xt] = fileparts(filename);
            fname = [fname '_horizon'];
            fn = fullfile(fpath,[fname xt]);
        end
        assert(~isempty(GeoKeyDirectoryTag),...
            'output is .tif so must include a ''GeoKeyDirectoryTag''')
        geotiffwrite(fn,horizons,R,'GeoKeyDirectoryTag',GeoKeyDirectoryTag)
        
        % same for distances
        if exist('distances','var')
            % change filename to include 'distance'
            [fpath,fname,xt] = fileparts(fn);
            if contains(fname,'horizon','IgnoreCase',true)
                fname = replace(lower(fname),'horizon','');
            elseif contains(lower(fname),'horz','IgnoreCase',true)
                fname = replace(lower(fname),'horz','');
            elseif contains(fname,'hor','IgnoreCase',true)
                fname = replace(lower(fname),'hor','');
            end
            fname = [fname '_distance'];
            fname = replace(fname,'__','_');
            fn = fullfile(fpath,[fname xt]);
            geotiffwrite(fn,distances,R,'GeoKeyDirectoryTag',GeoKeyDirectoryTag)
        end
    end

    function writeHDF5()
        % horizonsDF specific stuff
        h5file = filename;
        deflateLevel = 9; % compression factor
        group = '/Grid';
        horizonScale = 100;
        ihorizonsorz = float2integer(horizons,horizonScale,0,'int16');
        hfill = intmin('int16');
        
        % create file and write data
        dname = 'horizons';
        h5create(h5file,[group '/' dname],size(ihorizonsorz),...
            'Deflate',deflateLevel,...
            'ChunkSize',[size(ihorizonsorz,1) size(ihorizonsorz,2) 1],...
            'FillValue',hfill,...
            'DataType',class(ihorizonsorz));
        h5write(h5file,[ group '/' dname],ihorizonsorz)
        h5writeatt(h5file,[ group '/' dname],'Divisor',horizonScale);
        h5writeatt(h5file,[ group '/' dname],'Units','degree');
        if exist('distances','var')
            dname = 'distances';
            distanceScale = 1;
            if contains(class(distances),'single')
                dfill = single(NaN);
            else
                dfill = NaN;
            end
            h5create(h5file,[group '/' dname],size(distances),...
                'Deflate',deflateLevel,...
                'ChunkSize',[size(distances,1) size(distances,2) 1],...
                'FillValue',dfill,...
                'DataType',class(distances));
            h5write(h5file,[ group '/' dname],distances)
            h5writeatt(h5file,[ group '/' dname],'Divisor',distanceScale);
            h5writeatt(h5file,[ group '/' dname],'Units','meter');
        end
        
        % write Global Attributes
        if exist('proj','var')
            PS = whos('proj');
            h5create(h5file,['/' 'MapProjection'],PS.bytes)
            h5writeProjection(h5file,['/' 'MapProjection'],proj)
        end
        if exist('geoid','var')
            assert(isnumeric(geoid),...
                'geoid must be numeric pair [SemimajorAxis Eccentricity]')
            
            h5writeatt(h5file,'/','geoid',geoid)
        end
        h5writeatt(h5file,group,'azimuths',azm(:).');
        sz = GetSize(R);
        h5create(h5file,['/' 'SpatialRef'],sz)
        h5writeSpatialRef(h5file, ['/' 'SpatialRef'],R)
    end
end