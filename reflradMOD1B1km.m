function [ S ] = reflradMOD1B1km( file, varargin )
% [ S ] = reflradMOD1B1km( file [,bands] )
%extract refectances (for solar bands) and/or radiances (for emissive bands)
% and ancillary information from MODIS Level 1B 1KM HDF file
%
%Input
%   file - MODIS L1B HDF filename
%Optional input
%   bands - numeric vector (or cell vector of bands as characters) of band
%       numbers to extract
%
%Output
%   S - structure with reflectance (for solar bands) or radiance (for emissive
%       bands and ancillary information, all at same size of radiance images
%       (values are single-precision floats, with invalid values set to NaN)
%
%       Radiance grids have the name radiance_bandNN
%       Reflectance grids have the name reflectance_bandNN
%       Other values are:
%       Latitude (degrees)
%       Longitude (degrees)
%       SensorZenith (degrees)
%       SensorAzimuthS (degrees, 0 South, +ccw)
%       SolarZenith (degrees)
%       SolarAzimuthS (degrees, 0 South, +ccw)

% parse input
p = inputParser;
defaultBands = 1:36;
validationFcn = @(x) iscell(x) || ischar(x) || isnumeric(x);
addRequired(p,'file',@ischar)
addOptional(p,'bands',defaultBands,validationFcn)
parse(p,file,varargin{:});
file = p.Results.file;
if isnumeric(p.Results.bands)
    bands = p.Results.bands;
elseif iscell(p.Results.bands)
    bands = zeros(length(p.Results.bands),1);
    for k=1:length(bands)
        bands(k) = str2double(p.Results.bands{k});
    end
elseif ischar(p.Results.bands)
    bands = str2double(p.Results.bands);
end

% file name, without folder
S.file = fnstart(file);
assert(~isempty(strfind(S.file,'1KM')),[S.file ' not a 1KM file'])
x = hdfread(file,'Earth-Sun Distance');
S.EarthSunDistance = x{1};

% times
T = getTimes(file);
v = fieldnames(T);
for k=1:length(v)
    S.(v{k}) = T.(v{k});
end

% read the metadata values, but not the image data themselves
sds = {'EV_250_Aggr1km_RefSB','EV_500_Aggr1km_RefSB','EV_1KM_RefSB','EV_1KM_Emissive'};
ReflectanceScales = cell(1,length(sds)-1);
ReflectanceOffsets = cell(size(ReflectanceScales));
ValidRange = zeros(length(sds),2,'uint16');
allBands = cell(length(sds),1);
for d=1:length(sds)
    v = hdfread(file,[sds{d} '/valid_range']);
    ValidRange(d,:) = v{1};
    v = hdfread(file,[sds{d} '/band_names']);
    tb = textscan(replace(v{1}',',',' '),'%s');
    allBands{d} = tb{1};
    if d<=3
        ReflectanceScales{d} = hdfread(file,[sds{d} '/reflectance_scales']);
        ReflectanceOffsets{d} = hdfread(file,[sds{d} '/reflectance_offsets']);
    else
        % radiances only for the emissive bands
        RadianceScales = hdfread(file,[sds{d} '/radiance_scales']);
        RadianceOffsets = hdfread(file,[sds{d} '/radiance_offsets']);
    end
end

% read and convert just the bands we specify
holdSize = [];
for d=1:length(sds)
    % read this SDS if it overlaps with any of the bands specified
    readThisOne = false;
    theseBands = allBands{d}';
    for n=1:length(theseBands)
        for b=1:length(bands)
            if strncmp(theseBands{n},num2str(bands(b)),2)
                readThisOne = true;
                break;
            end
        end
        if readThisOne
            break;
        end
    end
    if readThisOne % at least one matching band
        U = hdfread(file,sds{d});
        for b=1:length(bands)
            for n=1:length(theseBands)
                if strncmp(theseBands{n},num2str(bands(b)),2)
                    if d<=3
                        Scale = ReflectanceScales{d}{1}(n);
                        Offset = ReflectanceOffsets{d}{1}(n);
                        var = ['reflectance_band' num2str(bands(b),'%02d')];
                    else
                        Scale = RadianceScales{1}(n);
                        Offset = RadianceOffsets{1}(n);
                        var = ['radiance_band' num2str(bands(b),'%02d')];
                    end
                    
                    raw = squeeze(U(n,:,:));
                    if isempty(holdSize)
                        holdSize = size(raw);
                    end
                    X = Scale*(single(raw)-Offset);
                    vx = ValidRange(d,:);
                    X(raw<min(vx)|raw>max(vx)) = NaN;
                    S.(var) = X;
                end
            end
        end
        % set reflectance bands empty if at night
        if strcmpi(S.DayNightFlag,'Night')
            fnames = fieldnames(S);
            for k=1:length(fnames)
                if strfind(fnames{k},'reflectance_band')
                    S.(fnames{k}) = [];
                end
            end
        end
    end
end

% ancillary information
A = ancillaryMOD1B1km(file,holdSize);
fn = fieldnames(A);
for k=1:length(fn)
    S.(fn{k}) = A.(fn{k});
end
end

% removes folder information, through final slash
function s=fnstart(c)
k = strfind(c,'\');
if isempty(k)
    k = strfind(c,'/');
end
if isempty(k)
    s = c;
else
    s = c(k(end)+1:end);
end
end

% ancillary data
function S = ancillaryMOD1B1km(file,imageSize)
% S = ancillaryMOD1B1km(file,imageSize)
%height, range, sensor and solar angles, and latitude and longitude from
%MODIS Level 1B 1Km file, interpolated to image of the radiance/reflectance
%values
%
%Input
%   file - MODIS filename
%
%Output (angles in degrees)
%   S - structure with elements in single precision:
%       SensorZenith,
%       SensorAzimuthS,
%       SolarZenith,
%       SolarAzimuthS,
%       Latitude,
%       Longitude,
%       Height (of land surface),
%       Range (sensor to pixel)
% Output azimuths are from SOUTH, positive ccw

loc = {'Latitude','Longitude'};
ang = {'SensorZenith','SensorAzimuth','SolarZenith','SolarAzimuth'};
dis = {'Height','Range'};

% latitudes and longitudes
for k=1:length(loc)
    X = hdfread(file,loc{k});
    v = hdfread(file,[loc{k} '/valid_range']);
    vR = v{1};
    X(X<min(vR) | X>max(vR)) = NaN;
    R.(loc{k}) = X;
end

% solar/view angles
for k=1:length(ang)
    X = hdfread(file,ang{k});
    v = hdfread(file,[ang{k} '/valid_range']);
    vR = v{1};
    v = hdfread(file,[ang{k} '/scale_factor']);
    sF = v{1};
    thisAng = single(X)*sF;
    thisAng(X<min(vR)|X>max(vR)) = NaN;
    if strfind(ang{k},'Zenith')
        R.(ang{k}) = thisAng;
    elseif strfind(ang{k},'Azimuth') % change to ccw from S
        if thisAng>=0
            R.(ang{k}) = 180-thisAng;
        else
            R.(ang{k}) = -(thisAng+180);
        end
    else
        error('ang %s not recognized',ang{k})
    end
end

% range and height
for k=1:length(dis)
    X = hdfread(file,dis{k});
    v = hdfread(file,[dis{k} '/valid_range']);
    vR = v{1};
    thisDis = single(X);
    thisDis(X<min(vR) | X>max(vR)) = NaN;
    R.(dis{k}) = thisDis;
end

% interpolate to the size of the image
fn = fieldnames(R);
dim = {'/line_numbers', '/frame_numbers'};
startVal = zeros(2,1);
diffVal = zeros(2,1);
[fullR,fullC] = ndgrid(1:imageSize(1),1:imageSize(2));
for k=1:length(fn)
    for d=1:length(dim)
        n = hdfread(file,[fn{k} dim{d}]);
        N = replace(n{1}',',',' ');
        kN = strfind(N,' ');
        startVal(d) = str2double(N(1:kN(1)));
        if length(kN)==1
            endVal = str2double(N(kN(1):end));
        else
            endVal = str2double(N(kN(1):kN(2)));
        end
        diffVal(d) = endVal-startVal(d);
    end
    % rows and columns in full image
    rows = startVal(1):diffVal(1):imageSize(1);
    cols = startVal(2):diffVal(2):imageSize(2);
    [r,c] = ndgrid(rows,cols);
    % use linear interpolation, but also use nearest-neighbor to keep NaNs
    % from propagating
    Fl = griddedInterpolant(r,c,R.(fn{k}));
    X = Fl(fullR,fullC);
    if any(isnan(X(:)))
        Fn = griddedInterpolant(r,c,R.(fn{k}),'nearest');
        Xn = Fn(fullR,fullC);
        t = isnan(X);
        X(t) = Xn(t);
    end
    S.(fn{k}) = X;
end
end


% times of Equator crossing and range
function T = getTimes(file)
I = hdfinfo(file);
% find core metadata attribute
for k=1:length(I.Attributes)
    if strncmpi(I.Attributes(k).Name,'CoreMetadata',12)
        A = I.Attributes(k);
        break;
    end
end
V = A.Value;
% variables
var = {'DayNightFlag','EquatorCrossingTime','RangeBeginningTime','RangeEndingTime'};
for k=1:length(var)
    if k==1
        ks = strfind(V,upper(var{k}));
        vs = V(ks(1):ks(2));
        kq = strfind(vs,'"');
        value = vs(kq(1)+1:kq(2)-1);
        T.(var{k}) = value;
    else % various times
        tvar = upper(var{k});
        dvar = upper(replace(var{k},'Time','Date'));
        kd = strfind(V,dvar);
        kt = strfind(V,tvar);
        vd = V(kd(1):kd(2));
        vt = V(kt(1):kt(2));
        % date
        kd = strfind(vd,'VALUE');
        vd = vd(kd(1):end);
        kq = strfind(vd,'"');
        valueD = vd(kq(1)+1:kq(2)-1);
        % time of day
        kt = strfind(vt,'VALUE');
        vt = vt(kt(1):end);
        kq = strfind(vt,'"');
        valueT = vt(kq(1)+1:kq(2)-1);
        t = datetime([valueD ' ' valueT],'InputFormat','yyyy-MM-dd HH:mm:ss.SSSSSS','TimeZone','UTC');
        t.Format = 'uuuu-MMM-dd HH:mm:ss ZZZZ';
        T.(var{k}) = t;
    end
end
end