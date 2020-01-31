function [ sunrise, sunset ] = sunRiseSet( lat,lon,matdates,varargin )
% [ sunrise, sunset ] = sunRiseSet( lat,lon,matdates [,P,T] )
%sunRiseSet UTC times of sunrise and sunset
%
% Input
%   lat,lon - latitudes and longitudes, must be same size if not scalars
%   matdates - dates in MATLAB datenum format, UTC but part of day ignored,
%       must be same size as lat,lon if not scalar and they are not scalar
% Optional input, must have both if used, and causes refraction to be
% considered
%   P,T - pressure (kPa) and temperature (K)

% arguments
narginchk(3,5)
p = inputParser;
addRequired(p,'lat',@isnumeric)
addRequired(p,'lon',@isnumeric)
addRequired(p,'matdates',@isnumeric)
addOptional(p,'P',[],@isnumeric)
addOptional(p,'T',[],@isnumeric)
parse(p,lat,lon,matdates,varargin{:})
P = p.Results.P;
T = p.Results.T;
assert(~xor(isempty(P),isempty(T)),...
    'if P or T are used, both must be specified')
assert(isequal(size(lat),size(lon)),...
    'if lat and lon not scalars, must be same size')
assert((isscalar(lat) || isscalar(matdates)) ||...
    isequal(size(lat),size(matdates)),...
    'if lat/lon and matdates not scalars, must be same size')

% make same size if one is a scalar
if xor(isscalar(lat),isscalar(matdates))
    if isscalar(lat)
        lat = ones(size(matdates))*lat;
        lon = ones(size(matdates))*lon;
    else
        matdates = ones(size(lat))*matdates;
    end
end
if ~isempty(P)
    assert(isscalar(P) || isequal(size(P),size(lat)),...
        'if P is not scalar, must be same size as lat/lon or matdates')
    if isscalar(P) && ~isscalar(lat)
        P = ones(size(lat))*P;
        T = ones(size(lat))*T;
    end
end
% hold size in case we vectorize
N = size(lat);
lat = lat(:);
lon = lon(:);
matdates = matdates(:);
refrac = ~isempty(P);
if refrac
    P = P(:);
    T = T(:);
end

r = zeros(size(lat));
s = zeros(size(r));
for k=1:length(lat)
    if refrac
        [r(k),s(k)] = sunUpDown(lat(k),lon(k),matdates(k),P(k),T(k));
    else
        [r(k),s(k)] = sunUpDown(lat(k),lon(k),matdates(k));
    end
end

% resize if necessary
if ~isequal(N,size(lat))
    sunrise = reshape(r,N);
    sunset = reshape(s,N);
else
    sunrise = r;
    sunset = s;
end

end

function [r,s] = sunUpDown(lat,lon,d,varargin)
if nargin>3
    P = varargin{1};
    T = varargin{2};
    refrac = true;
else
    refrac = false;
end

[declin,~,~,eqT] = Ephemeris(d);
solarNoon = (720-4*lon-eqT)/1440;
arg = cosd(90.833)./(cosd(lat).*cosd(declin))-tand(lat).*tand(declin);
arg(arg<-1) = -1;
arg(arg>1) = 1;
HA = acosd(arg); % hour angle, sunrise
r = solarNoon-HA*4/1440;
s = solarNoon+HA*4/1440;
end
