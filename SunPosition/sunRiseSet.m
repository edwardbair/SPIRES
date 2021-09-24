function [ sunrise, sunset ] = sunRiseSet( lat,lon,matdates,varargin )
% [ sunrise, sunset ] = sunRiseSet( lat,lon,matdates [,P,T] )
%sunRiseSet times of sunrise and sunset, as MATLAB datetimes
%
% Input
%   lat,lon - latitudes and longitudes, must be same size if not scalars
%   matdates - dates in MATLAB datenum or datetime format, hour and
%       minutes are ignored
%       if datetime format, specify the TimeZone, otherwise UTC is assumed
%       must be scalar or same size as lat,lon
% Optional input, must have both if used, and causes refraction that depends
%   on pressure and temperature to be considered
%   P,T - pressure (hPa) and temperature (K)
%
% Output
%   MATLAB datetimes corresponding to the input dates in the time zone
%       specified in the input
%   if 24 hr daylight, sunrise and sunset are set to the beginning and
%       ending of the day
%   if 24 hr darkness, sunrise and sunset are both set to NaT

% arguments
narginchk(3,5)
p = inputParser;
addRequired(p,'lat',@(x) isnumeric(x) && all(abs(x(:))<=90))
addRequired(p,'lon',@(x) isnumeric(x) && all(abs(x(:))<=180))
addRequired(p,'matdates',@(x) isnumeric(x) || isdatetime(x))
addOptional(p,'P',[],@isnumeric)
addOptional(p,'T',[],@isnumeric)
parse(p,lat,lon,matdates,varargin{:})
P = p.Results.P;
T = p.Results.T;
assert(~xor(isempty(P),isempty(T)),...
    'if P or T are used, both must be specified')
lat = p.Results.lat;
lon = p.Results.lon;
if isdatetime(p.Results.matdates)
    matdates = p.Results.matdates;
    if isempty(matdates.TimeZone)
        matdates.TimeZone = 'Etc/GMT';
    end
    % eliminate the part of the day
    dn = datenum(matdates);
    matdates = datetime(floor(dn),'ConvertFrom','datenum',...
        'TimeZone',matdates.TimeZone);
else
    matdates = datetime(floor(p.Results.matdates),'ConvertFrom','datenum',...
        'timezone','Etc/GMT');
end
assert((isscalar(lat) || isscalar(matdates)) ||...
    isequal(size(lat),size(matdates)),...
    'if lat/lon and matdates not scalars, must be same size')

% make same size if one is a scalar
[lat,lon,matdates] = checkSizes(lat,lon,matdates);
if ~isempty(P)
    assert(isscalar(P) || isequal(size(P),size(lat)),...
        'if P is not scalar, must be same size as lat/lon or matdates')
end
% hold size in case we vectorize
N = size(lat);
lat = lat(:);
lon = lon(:);
matdates = matdates(:);

% values just based on declination at time zone at noon (not accounting for
% change in declination during the day)
[declin,~,~,eqT] = EarthEphemeris(matdates+hours(12));
% refraction correction
refrac = ~isempty(P);
if refrac
    P = P(:);
    T = T(:);
end
if refrac
    if isscalar(P)
        muR = fzero(@xref,0);
    else
        assert(isequal(size(P),size(lat)),'P and lat must be of same sizes')
        muR = zeros(size(lat));
        % P & T are passed to xref internally, so must do one by one
        holdP = P;
        holdT = T;
        for k=1:length(holdP)
            P = holdP(k);
            T = holdT(k);
            muR(k) = fzero(@xref,0);
        end
    end
else
    muR = 0;
end
arg = muR./(cosd(lat).*cosd(declin))-tand(lat).*tand(declin);
r = NaT(size(lat),'TimeZone',matdates.TimeZone);
s = NaT(size(lat),'TimeZone',matdates.TimeZone);
ts = arg<-1; %total sunlight
td = arg>1; %total darkness, leave as NaT
t = ~ts & ~td; %all else
tz = hoursFromTimeZone(matdates);
solarNoon = (720-4*lon-eqT+tz*60)/1440;
if nnz(t)
    HA = acosd(arg(t));
    r(t) = matdates(t)+days(solarNoon(t)-4*HA/1440);
    s(t) = matdates(t)+days(solarNoon(t)+4*HA/1440);
end
if nnz(ts)
    HA = 180;
    r(ts) = matdates(ts)+days(solarNoon(ts)-4*HA/1440);
    s(ts) = r(ts)+days(1);
end

% resize if necessary
if ~isequal(N,size(lat))
    sunrise = reshape(r,N);
    sunset = reshape(s,N);
else
    sunrise = r;
    sunset = s;
end
    function mu_Refracted = xref(mu0)
        mu_Refracted = refracted(mu0,P,T);
    end
end

function hz=hoursFromTimeZone(d)
assert(isdatetime(d),'argument must be a datetime')
if isempty(d.TimeZone)
    hz = 0;
else
    d2 = d;
    d2.TimeZone = 'Etc/UTC';
    [yy,mm,dd,hh,m,s] = datevec(d);
    tempD = datetime(yy,mm,dd,hh,m,s,'TimeZone',d2.TimeZone);
    [H,M,~] = hms(tempD-d2);
    hz = H+M/60;
end
end