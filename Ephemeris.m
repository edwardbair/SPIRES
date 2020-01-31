function [ declin, radiusvector, solar_lon, varargout ] = Ephemeris( matdates )
%Ephemeris: [ declin, radiusvector, solar_lon [,eqTime] ] = Ephemeris( matdates )
% using calculations from the NOAA Solar Calculator
%   http://www.esrl.noaa.gov/gmd/grad/solcalc/
%
%input
%   matdates - scalar or vector (or matrix) of MATLAB datenums (UTC)
%
%output
%   declination & solar longitude in degrees
%   radius vector in AU
%   optional - equation of time

% matdates to julian centuries
assert(isnumeric(matdates),'argument must be numeric')
dt = datetime(matdates(:),'ConvertFrom','datenum');
G = (juliandate(dt)-2451545)/36525; % Julian century

% calculations - the letters correspond to the collumns in the spreadsheets
% at the NOAA Solar Calculator
I = mod(280.46646+G.*(36000.76983 + G*0.0003032),360); % geom mean long Sun, degrees
J = 357.52911+G.*(35999.05029 - 0.0001537*G); % geom mean Sun anomaly, degrees
K = 0.016708634-G.*(0.000042037+0.0000001267*G); % Earth orbit eccentricity
L = sind(J).*(1.914602-G.*(0.004817+0.000014*G))+sind(2*J).*...
    (0.019993-0.000101*G)+sind(3*J)*0.000289; % Sun eq of Ctr
M = I+L; % Sun true longitude
N = J+L; % Sun true anomaly
radiusvector = (1.000001018*(1-K.^2))./(1+K.*cosd(N));
P = M-0.00569-0.00478*sind(125.04-1934.136*G); % Sun apparent long
Q = 23+(26+((21.448-G.*(46.815+G.*(0.00059-G*0.001813))))/60)/60; % mean obliq ecliptic
R = Q+0.00256*cosd(125.04-1934.136*G); % obliquity correction
T = asind(sind(R).*sind(P)); % declination
U = tand(R/2).^2; % var y
V = 4*rad2deg(U.*sind(2*I)-2*K.*sind(J)+4*K.*U.*sind(J).*cosd(2*I)-...
    0.5*U.^2.*sind(4*I)-1.25*K.^2.*sind(2*J)); % eq of time
E = matdates(:)-floor(matdates(:)); % local time
AB = mod(E*1440+V,1440); % true solar time, min
t = AB<0;
solar_lon = AB/4;
solar_lon(t) = -(solar_lon(t)+180);
solar_lon(~t) = -(solar_lon(~t)-180);
declin = T;

% revert to matrix, if needed
if ~isequal(size(matdates),size(declin))
    declin = reshape(declin,size(matdates));
    solar_lon = reshape(solar_lon,size(matdates));
    radiusvector = reshape(radiusvector,size(matdates));
    V = reshape(V,size(matdates));
end

% equation of time
if nargout>0
    varargout{1} = V;
end

end