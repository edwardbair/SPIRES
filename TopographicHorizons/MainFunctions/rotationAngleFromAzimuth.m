function [angToRotate,fwdAz,backAz] = rotationAngleFromAzimuth(desiredAzimuth,rasterref,varargin)
% [angToRotate,fwdAz,backAz] = rotationAngleFromAzimuth(desiredAzimuth,rasterref,varargin)
%estimate angle to rotate grid to get desired azimuth
%returned values provide both forward and backward azimuths (because they
%can be calculated on the same pass through horizonAlongProfile)
%
%Input
% desiredAzimuth - degrees, either +/-180 or 0 to 360, depending on how
%   function azimuthPreference is set
% rasterref - of the input grid, can be either Map or Geographic
%Third input, name/value pair
% if rasterref is a Map object
%   'projection' - structure describing the projection to convert from map to
%       geographic coordinates
% if rasterref is a Geographic object
%   'planet' - default 'wgs84'
%
%Output
% angToRotate - angle to rotate the grid to get the desired azimuth,
%   restricted to the -90 to +90 range
% fwdAz and backAz - forward and backward azimuths, degrees, either +/-180
%   or 0 to 360 depending on return value from function azimuthPreference

p = inputParser;
addRequired(p,'desiredAzimuth',@(x) isscalar(x) && isnumeric(x) &&...
    x>=-180 && x<=360)
addRequired(p,'rasterref',@(x) contains(class(x),'rasterref'))
addParameter(p,validatestring('projection',{'proj','projection'}),...
    struct([]),@(x) isstruct(x) || isempty(x))
addParameter(p,'planet','wgs84',@ischar)
parse(p,desiredAzimuth,rasterref,varargin{:})

isGeographic = contains(class(rasterref),'geographic','IgnoreCase',true);
if isGeographic
    if ~isempty(rasterref.GeographicCRS)
        E = rasterref.GeographicCRS.Spheroid;
    else
        E = referenceEllipsoid(p.Results.planet);
    end
else
    if ~isempty(rasterref.ProjectedCRS)
        proj = rasterref.ProjectedCRS;
        E = rasterref.ProjectedCRS.GeographicCRS.Spheroid;
    else
        proj = p.Results.projection;
        if isstruct(proj)
            E = oblateSpheroid;
            E.SemimajorAxis = proj.geoid(1);
            E.Eccentricity = proj.geoid(2);
        else
            E = referenceEllipsoid(p.Results.planet);
        end
    end
    assert(~isempty(proj),'''projection'' is needed because rasterref is a Map object without a ProjectedCRS')
end

% because MATLAB distance function returns azimuths in 0 to 360 range, convert
if azimuthPreference
    azm = 180-desiredAzimuth;
else
    azm = desiredAzimuth;
end
% restrict to range 90 to 270, as backward calculation will get the other direction
if azm<90
    azm = azm+180;
elseif azm>=270
    azm = azm-180;
end

% create lat-lon grid, just 10% in size (but don't do for small grids, less
% than 500 pixels)
R = rasterref;
if min(R.RasterSize)>500
    gridSize = round(R.RasterSize/10);
else
    gridSize = R.RasterSize;
end
if isGeographic
    [lon,lat] = meshgrid(linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),gridSize(2)),...
        linspace(R.LatitudeLimits(2),R.LatitudeLimits(1),gridSize(1)));
else
    [x,y] = meshgrid(linspace(R.XWorldLimits(1),R.XWorldLimits(2),gridSize(2)),...
        linspace(R.YWorldLimits(2),R.YWorldLimits(1),gridSize(1)));
    [lat,lon] = projinv(proj,x,y);
end

% get good start guess
rot0 = guessRotation(lat,lon,azm,E);

% solve for rotation angle
[angToRotate,solvedDiff,exitflag,output] = fzero(@azmDiff,rot0); %#ok<ASGLU>
assert(exitflag==1,'fzero did not converge on a solution')
solvedAz = solvedDiff+azm;
if azimuthPreference
    solvedAz = 180-solvedAz;
end
fwdAz = solvedAz;
backAz = solvedAz-180;

    function rotDiff=azmDiff(rotAng)
        azmRot = meanAzimuth(rotAng);
        rotDiff = (azmRot-azm);
    end

    function azmR = meanAzimuth(rotationAngle)
        az = [];
        latR = lat;
        lonR = lon;
        latR(latR==0) = -Inf;
        lonR(lonR==0) = -Inf;
        latR = imrotate(latR,rotationAngle);
        lonR = imrotate(lonR,rotationAngle);
        lonR(lonR==0) = NaN;
        latR(latR==0) = NaN;
        lonR(isinf(lonR)) = 0;
        latR(isinf(latR)) = 0;
        mididx = round(size(lonR,2)/2);
        mididx = [mididx-20 mididx+20];
        for m=1:length(mididx)
            lonProfile = lonR(:,mididx(m));
            latProfile = latR(:,mididx(m));
            lonProfile = lonProfile(~isnan(lonProfile));
            latProfile = latProfile(~isnan(latProfile));
            [D,azhold] = distance(latProfile(1),lonProfile(1),latProfile,lonProfile,E); %#ok<ASGLU>
            az = cat(1,az,azhold);
        end
        azmR = atan2d(median(sind(az),'omitnan'),median(cosd(az),'omitnan'));
        if azmR<0
            azmR = 360+azmR;
        end
    end
end

function rot0=guessRotation(lat,lon,azm,E)
tryRotation = -90:30:90;
% imrotate fills with zeros, so fill zero lats and lons with -Inf
lat(lat==0) = -Inf;
lon(lon==0) = -Inf;
tryAzm = zeros(size(tryRotation));
for n=1:length(tryRotation)
    latR = imrotate(lat,tryRotation(n));
    lonR = imrotate(lon,tryRotation(n));
    latR(latR==0) = NaN;
    lonR(lonR==0) = NaN;
    latR(isinf(latR)) = 0;
    lonR(isinf(lonR)) = 0;
    % distance and azimuth in middle 11 longitudes
    mididx = round(size(lonR,2)/2);
    mididx = [mididx-5 mididx+5];
    xAzm = [];
    for m=1:length(mididx)
        idx = mididx(m);
        if mididx(m)<1
            idx = 1;
        end
        if mididx(m)>size(latR,2)
            idx = size(latR,2);
        end
        latProfile = latR(:,idx);
        if mididx(m)>size(lonR,2)
            idx = size(lonR,2);
        end
        lonProfile = lonR(:,idx);
        latProfile = latProfile(~isnan(latProfile));
        lonProfile = lonProfile(~isnan(lonProfile));
        [D,az] = distance(latProfile(1),lonProfile(1),latProfile,lonProfile,E); %#ok<ASGLU>
        xAzm = cat(1,xAzm,az);
    end
    tryAzm(n) = atan2d(median(sind(xAzm),'omitnan'),median(cosd(xAzm),'omitnan'));
end
% find where mean azimuth changes sign
tz = tryAzm<0;
tryAzm(tz) = 360+tryAzm(tz);
tryAzm = tryAzm-azm;
rot0 = [];
for n=2:length(tryAzm)
    if tryAzm(n-1)*tryAzm(n)<=0
        rot0 = [tryRotation(n-1) tryRotation(n)];
        break
    end
end
assert(~isempty(rot0),'cannot find rotation angles spanning solution')
end