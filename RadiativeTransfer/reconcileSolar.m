function [ estSolar,mu0,airmass ] = reconcileSolar(inputSolar,meastime,zoneOffset,latitude,longitude,varargin )
% [ estSolar ] = reconcileSolar(inputSolar,meastime,zoneOffset,latitude,longitude,... )
%make sure values of solar radiation are reasonable, based on location and
%time of day, otherwise correct from estimated optical depth
%
% Input
%   inputSolar - measured solar radiation on flat surface (scalar or
%       vector)
%   meastime - time of measurement in MATLAB daynum format (scalar or
%       vector, same size as inputSolar)
%   zoneOffset - time zone offset in hours (e.g., PST is +8 hr)
%   latitude - location of measurement (degrees)
%   longitude (degrees)
% Optional input, name-value pairs
%   'S0' - solar constant, W/m^2 default 1368
%   'T' - air temperature, deg K (default 288)
%   'P' - air pressure, kPa (default 100)
%   'Z' - elevation, m (default 0)
%   (T & P can be vectors of same size as inputSolar)
%   (if Z is specified and P is not, scale P hydrostatically)
%
% Output
%   estSolar - estimate correct value of solar radiation
%   mu0 - cosine of solar angle
%   airmass - by Kasten's equation

p = inputParser;
defaultS0 = 1368;
defaultP = 100;
defaultT = 288;
defaultZ = 0;

addRequired(p,'inputSolar',@isnumeric);
addRequired(p,'meastime',@isnumeric);
addRequired(p,'zoneOffset',@isnumeric);
addRequired(p,'latitude',@isnumeric);
addRequired(p,'longitude',@isnumeric);
addParameter(p,'S0',defaultS0,@isnumeric);
addParameter(p,'P',defaultP,@isnumeric);
addParameter(p,'T',defaultT,@isnumeric);
addParameter(p,'Z',defaultZ,@isnumeric);

parse(p,inputSolar,meastime,zoneOffset,latitude,longitude,varargin{:});

% assign variables if used more than once (otherwise, use p structure)
S0 = p.Results.S0;
P = p.Results.P;
TK = p.Results.T;
Z = p.Results.Z;

[declin,radiusvector,solar_lon] = Ephemeris(meastime+p.Results.zoneOffset/24);

% adjust pressure if elevation specified but pressure not
if ~isequal(Z,defaultZ) && isequal(P,defaultP)
    P = AirPressure(defaultP,defaultT,defaultZ,TK,Z);
end
[mu0,~,airmass] = sunang(p.Results.latitude,p.Results.longitude,declin,...
    solar_lon,true,P,TK);

% optical depths, based on measurements
threshold = 0.25;
tau = (log(mu0*S0./(radiusvector.^2.*inputSolar)))./airmass;
tau(mu0<threshold) = NaN;
badTau = (tau<=0 | isnan(tau) | isinf(tau) | ~isreal(tau)) & mu0>0;
tau(badTau) = NaN;
tnan = isnan(tau);
mintau = prctile(tau(~tnan),.1);
maxtau = prctile(tau(~tnan),99.9);
tokay = ~tnan & tau>=mintau & tau<=maxtau;
F = fit(meastime(tokay),tau(tokay),'smoothingspline','Weights',mu0(tokay));
that = truncateLimits(F(meastime),mintau,maxtau);
that(mu0<=0) = NaN;

estSolar = S0.*mu0.*radiusvector.^(-2).*exp(-that.*airmass);
estSolar(isnan(that) | isnan(estSolar) | estSolar<0) = 0;

end