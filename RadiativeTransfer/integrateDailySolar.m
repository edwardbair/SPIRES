function [ value,timeInterval ] = integrateDailySolar( lat,lon,matdate,timezone,varargin )
% [ value,timeInterval ] = integrateDailySolar( lat,lon,matdate,timezone,... )
%integrate solar radiation over a day
%
% Input
%   lat,lon - of location, scalar
%   matdate - date for calculation, scalar
%   timezone - hours, scalar, - to West of UTC, + to East
% Optional input, name-value pairs, case insensitive and characters beyond
%   the 4th are ignored
%   'S0' - direct radiation normal to the beam, default 1.0 for scaling
%       (scalar)
%   'albedo' - in this case, the value returned is net radiation (scalar)
%   'elevation' - meters (scalar), default 0.0;
%   'slope','aspect' - slopes and aspects of locations, degrees, both must
%       be specified (aspects 0 deg to S, + to E, - to W), default 0.0
%   'horizon' - horizon vector, values being the sine of the angles to the
%       horizon, of length nHorz, where nHorz horizons cover the
%       directions from -180 (north) to 0 (south) to +180 (north), with
%       first and last values repeated (default no horizon)
%   'viewfactor' - sky view factor, default 1.0
%   For the following two sets of variables, if measurements are used, the
%       atmospheric properties are ignored
%   (following three vectors if using surface measurements to drive)
%       'meastime' - MATLAB datenums of the measurement times
%       'direct' - direct solar radiation values at those times
%       'diffuse' - diffuse solar radiation values
%   (following three vectors if using atmospheric model to drive, scalar or
%       same size as lat/lon)
%       'tau' - optical depth of atmosphere, default 0.0 in which case the other
%           two variables are ignored
%       'omega' - single scattering albedo, default 0.9
%       'asymmetry' - scattering asymmetry factor, default 0.8
%
% Output
%   value - incoming solar radiation, or net solar if albedo is specified
%   timeInterval - integration time start and finish, for local time zone

% parse input
narginchk(4,24)
[lat,lon,matdate,timezone,useHorz,useMeas,S0,albedo,elev,slope,aspect,...
    horz,viewF,meastime,direct,diffuse,tau,omega,g] =...
    parseInput(lat,lon,matdate,timezone,varargin{:});

% if measurements are used, fit an interpolant
if useMeas
    Fdir = fit(meastime,direct,'linearinterp');
    Fdif = fit(meastime,diffuse,'linearinterp');
else
    Fdir = [];
    Fdif = [];
end

% sunrise and set times for each, based on standard refraction, which will
% ensure that start and end times have zero solar radiation
[sunrise,sunset] = sunRiseSet(lat,lon,matdate);
sunrise = sunrise+matdate;
sunset = sunset+matdate;

% air pressure and temperature based on standard atmosphere (used to
% calculate refraction and airmass, so approximation is good enough)
P = AirPressure(100,288,0,[],elev);
lapse = -0.0065;
T = 288+elev*lapse;

% if horizons are used, fit an interpolant
if useHorz
    phi = linspace(-180,180,length(horz))';
    Fhorz = fit(phi,horz,'linearinterp');
else
    Fhorz = [];
end

if isempty(albedo)
    thisAlbedo = [];
else
    thisAlbedo = albedo;
end

if isnan(sunrise) && isnan(sunset)
    value = 0;
    timeInterval = nan(1,2);
else
    value = dailySolar(lat,lon,sunrise,sunset,S0,thisAlbedo,...
        P,T,slope,aspect,Fhorz,viewF,Fdir,Fdif,tau,omega,g);
    timeInterval = [sunrise sunset]+timezone/24;
end
end

function flag = visible(mu0,phi0,Fhorz)
% is this point visible?
assert(isequal(size(mu0),size(phi0)),'sizes of mu0 and phi0 must be same')
if isempty(Fhorz)
    flag = true(size(mu0));
else
    flag = mu0>Fhorz(phi0)';
end
end

function result=dailySolar(latitude,longitude,startTime,endTime,S0,albedo,...
    pressure,temperature,slope,aspect,Fhorz,viewF,Fdir,Fdif,tau,omega,g)
% integrates daily solar radiation for a specific location

% analytic solution if no slope and no atmosphere
% if slope==0 && tau==0
%     [declin,radiusvector,~] = Ephemeris(((startTime+endTime)/2)/86400);
%     result = 2*(acos(-tand(latitude)*tand(declin))*sind(latitude)*sind(declin)+...
%         cosd(latitude)*cosd(declin)*...
%         sqrt((1-tand(latitude)*tand(declin))*(1+tand(latitude)*tand(declin))));
%     result = result*S0/radiusvector^2;
%     return;
% end

result = integral(@instantSolar,startTime*86400,endTime*86400,...
    'AbsTol',1.e-8,'RelTol',1e-4);
if ~isempty(albedo)
    result = result*(1-albedo);
end

    function value=instantSolar(thisTime)
        thisTime = thisTime/86400;
        [declin,radiusvector,solarlon] = Ephemeris(thisTime);
        % sun angle and visibility
        if tau==0 && isempty(Fdir)
            [mu0, phi0, ~] = sunang(latitude,longitude,declin,solarlon);
        else
            [mu0, phi0, airmass] = sunang(latitude,longitude,declin,...
                solarlon,true,pressure,temperature);
        end
        if slope==0
            mu = mu0;
        else
            mu = sunslope(mu0,phi0,slope,aspect);
            mu(isnan(mu) | mu<0) = 0;
        end
        isVis = visible(mu0,phi0,Fhorz);
        if tau==0 && isempty(Fdir)
            value = (S0./radiusvector.^2).*mu;
            value(~isVis) = 0;
        elseif ~isempty(Fdir)
            diffuse = Fdif(thisTime);
            direct = Fdir(thisTime)/mu0;
            direct(mu0==0) = 0;
            value = direct*mu+diffuse*viewF;
        else
            assert(tau>0,'shouldn''t get to this statement if tau is zero')
            Trans = zeros(size(mu0));
            Beam = zeros(size(mu0));
            for k=1:length(mu0)
                if isnan(airmass(k))
                    Trans(k) = 0;
                    Beam(k) = 0;
                else
                    [~, Trans(k), Beam(k), ~] = twostream(mu0(k),omega,g,...
                        tau,'airmass',airmass(k)*pressure/100);
                end
            end
            direct = (S0./radiusvector.^2).*Beam;
            diffuse = (S0./radiusvector.^2).*Trans-direct;
            direct(~isVis) = 0;
            value = direct.*mu+diffuse.*viewF;
            scatter(thisTime,direct)
            hold on;
            scatter(thisTime,diffuse)
            scatter(thisTime,value)
        end
    end
end

function [lat,lon,matdate,timezone,useHorz,useMeas,S0,albedo,...
    elev,slope,aspect,horz,viewF,meastime,direct,diffuse,tau,omega,g] =...
    parseInput(LAT,LON,MATDATE,TIMEZONE,varargin)

% varargin switch to lower case and truncate to 4 chars
if ~isempty(varargin)
    for k=1:2:length(varargin)
        varargin{k} = lower(varargin{k});
        if length(varargin{k})>4
            s = varargin{k};
            varargin{k} = s(1:4);
        end
    end
end

p = inputParser;
vSfcn = @(x) validateattributes(x,{'numeric'},{'scalar','>=',-12,'<=',12});
vDfcn = @(x) isnumeric(x) && isscalar(x);
addRequired(p,'LAT',@isnumeric)
addRequired(p,'LON',@isnumeric)
addRequired(p,'MATDATE',vDfcn)
addRequired(p,'TIMEZONE',vSfcn)
% defaults
defaultS0 = 1;
defaultAlbedo = [];
defaultElevation = 0;
defaultSlope = 0;
defaultAspect = 0;
defaultHorizon = [];
defaultViewFactor = 1;
defaultMeastime = [];
defaultDirect = [];
defaultDiffuse = [];
defaultTau = 0;
defaultOmega = 0.9;
defaultAsymmetry = 0.8;
% optional arguments
addParameter(p,'s0',defaultS0,vDfcn)
addParameter(p,'albe',defaultAlbedo,vDfcn)
addParameter(p,'elev',defaultElevation,vDfcn)
addParameter(p,'slop',defaultSlope,vDfcn)
addParameter(p,'aspe',defaultAspect,vDfcn)
addParameter(p,'hori',defaultHorizon,@isnumeric)
addParameter(p,'view',defaultViewFactor,vDfcn)
addParameter(p,'meas',defaultMeastime,@isnumeric)
addParameter(p,'dire',defaultDirect,@isnumeric)
addParameter(p,'diff',defaultDiffuse,@isnumeric)
addParameter(p,'tau',defaultTau,vDfcn)
addParameter(p,'omeg',defaultOmega,vDfcn)
addParameter(p,'asym',defaultAsymmetry,vDfcn)
parse(p,LAT,LON,MATDATE,TIMEZONE,varargin{:});

% check sizes
lat = p.Results.LAT;
lon = p.Results.LON;
matdate = p.Results.MATDATE;
timezone = p.Results.TIMEZONE;

% optional arguments
S0 = p.Results.s0;
albedo = p.Results.albe;
elev = double(p.Results.elev);
slope = double(p.Results.slop);
aspect = double(p.Results.aspe);
viewF = double(p.Results.view);
tau = p.Results.tau;
omega = p.Results.omeg;
g = p.Results.asym;
meastime = p.Results.meas;
useMeas = ~isempty(meastime);
direct = p.Results.dire;
diffuse = p.Results.diff;
if useMeas
    assert(isequal(size(meastime),size(direct),size(diffuse)),...
        'if measurements specified, vectors of meastime, direct, and diffuse must be same size')
    % make sure measurement times cover the matdate
    mstart = floor(meastime);
    assert(~isempty(find(mstart==meastime,1)),...
        ['no measurement times for date ' num2str(meastime) ' ' datestr(meastime)])
end
horz = p.Results.hori;
useHorz = ~isempty(horz);
end