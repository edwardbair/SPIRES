function [horzAng,horzDis,horzAzm,varargout] = horizonAlongProfile(lat,lon,elev,refEllipsoid)
% [horzAng,horzDis,horzAzm [,horzPts]] = horizonAlongProfile(lat,lon,elev,refEllipsoid)
%
%calculates points that form the horizon points in the forward direction
%along an elevation profile
%(to go in the backward direction, just flip the input vectors)
%
% based on Dozier, J., J. Bruno, and P. Downey (1981), A faster solution to
% the horizon problem, Computers and Geosciences, 7, 145-151,
% doi: 10.1016/0098-3004(81)90026-1
%
% adapted here to lat-lon coordinates and accounting for planet size as
% specified in the reference ellipsoid
%
% INPUT - vectors of same length, ordered but need not be equispaced
%   lat - latitude vector of points
%   lon - longitude vector
%   elev - elevations
%   refEllipsoid - output from referenceEllipsoid or from
%       GeographicCRS.Spheroid in raster reference
%
% OUTPUT
%   horzAng - elevation angles to the horizons
%   horzDis - distance to horizon for each point (for points that are their
%       own horizons, distance to the end of the profile)
%   horzAzm - azimuth angles to the horizons
% Optional output
%   horzPts - indices of points that form the horizons

p = inputParser;
addRequired(p,'lat',@isvector)
addRequired(p,'lon',@isvector)
addRequired(p,'elev',@isvector)
addRequired(p,'refEllipsoid',@(x) contains(class(x),'referenceEllipsoid') ||...
    contains(class(x),'oblateSpheroid'))
parse(p,lat,lon,elev,refEllipsoid);

% check input
assert(isequal(size(lat),size(lon),size(elev)),...
    'input lat, lon, and elev must be vectors of same size');

t = ~isnan(elev);
if nnz(t)
    [horzPts,horzAng,horzDis,horzAzm] = horz1fPts(lat(t),lon(t),...
        double(elev(t)),refEllipsoid);
    % only points that are their own horizons have horizon angle zero
    % (set others with angle zero to NaN)
    tsame = horzPts==(1:length(horzPts)).';
    if any(horzAng(~tsame)==0)
        horzAng(~tsame & horzAng==0) = NaN;
    end
    if any(~t)
        % fill NaNs for all variables where elevation is NaN
        ha = nan(size(lat));
        hd = nan(size(ha));
        haz = nan(size(ha));
        hp = nan(size(ha));
        ha(t) = horzAng;
        hd(t) = horzDis;
        haz(t) = horzAzm;
        hp(t) = horzPts;
        horzAng = ha;
        horzDis = hd;
        horzAzm = haz;
        horzPts = hp;
    end
else
    horzPts = zeros(size(lat));
    horzAng = nan(size(lat));
    horzDis = nan(size(horzAng));
    horzAzm = nan(size(horzAng));
end
if nargout>3
    varargout{1} = horzPts;
end
end

function [Hpt,horzAng,horzDis,horzAzm] = horz1fPts(lat,lon,Z,E)
% 1D horizon function in the forward direction
% Hpt are indices of horizon points

% check size
if length(Z)<2
    horzAng = zeros(size(Z));
    horzDis = zeros(size(Z));
    horzAzm = nan(size(Z));
    Hpt = (1:length(Z)).';
    return
end

% distances along the profile based on lat-lon
[D,az] = distance(lat(1),lon(1),lat,lon,E);
% derivative at each point (to calculate slope when points are adjacent)
gradAlong = gradient(Z)./gradient(D);
gradAlong(gradAlong<0) = 0;

Hpt = zeros(size(Z));
assert(iscolumn(Hpt),'lat, lon, elev must be column vectors')
% end point is its own horizon
N = length(Hpt);
Hpt(end) = N;
holdSlope = zeros(size(Hpt));
holdDistance = zeros(size(Hpt));

% first azimuth is always zero (to itself) so fix
az(1) = az(2);
% to get median azimuth, allow for discontinuity across north, and skip
% first few as they tend to be noisy
if length(az)>=10
    az(1:4) = [];
end
azm = atan2d(nanmedian(sind(az)),nanmedian(cosd(az)));
% atan2d returns in range -180 to +180 whereas distance returns in range 0 to 360
if azm<0
    azm = azm+360;
end

% main loop
for ii=N-1:-1:1 % Loop from next-to-end backward to beginning
    % until slope from i to j is greater than or equal to the slope from
    % j to j's horizon (code uses ii and jj to not confuse with 1i or 1j).
    
    % Start with adjacent point
    jj = ii+1;
    % inner loop until horizon found
    found = false;
    while ~found
        sij = thisSlope(ii,jj); % slope i to j
        sjh = thisSlope(jj,Hpt(jj)); % slope j to its horizon
        if sij<sjh
            % try j's horizon next
            jj = Hpt(jj);
        else
            found = true;
            if sij>sjh % slope i to j greater than j to H(j) so j is i's horizon
                Hpt(ii) = jj;
                holdDistance(ii) = D(jj)-D(ii);
            elseif sij==0 % i is its own horizon
                Hpt(ii) = ii;
                holdDistance(ii) = D(end)-D(ii);
            else % slopes are equal so H(j) is i's horizon
                Hpt(ii) = Hpt(jj);
                holdDistance(ii) = D(Hpt(jj))-D(ii);
            end
        end
    end
    holdSlope(ii) = sij;
end
% convert slopes to angles in degrees
horzAng = atand(holdSlope);
% threshold distance and horizon angle at which we worry about planet
% curvature (100 km and 3 degrees)
horzDistanceThreshold = 100e3;
horzAngleThreshold = 3;
tdz = holdDistance>horzDistanceThreshold & horzAng>horzAngleThreshold;
if nnz(tdz)==0
    horzDis = holdDistance;
else
    [horzAng,horzDis,azm] =...
        invokeGeodetic(Hpt,horzAng,holdDistance,azm,lat,lon,Z,horzDistanceThreshold,E);
    if ~isscalar(azm)
        azm = atan2d(nanmedian(sind(azm)),nanmedian(cosd(azm)));
    end
    if azm<0
        azm = 360+azm;
    end
end
% azimuths either +/-180 or 0 to 360, see the azimuthPreference function
% (in the MATLAB distance function, output azimuths are 0 to 360)
if azimuthPreference
    horzAzm = 180-azm;
else
    horzAzm = azm;
end

    function s=thisSlope(kk,nn) % Z, D, gradAlong passed internally
        if Z(kk)>=Z(nn)
            s = 0;
        elseif nn-kk==1 % use central difference
            s = gradAlong(kk);
        else
            s = (Z(nn)-Z(kk))/(D(nn)-D(kk)); % forward difference
        end
    end
end

function [horzAng,horzDis,horzAzm] = invokeGeodetic(Hpt,ang,dis,azm,lat,lon,Z,threshold,E)
%returns geodetic distance, angles, and azimuths only for distance to horizon
%above threshold
%
% The code takes care of the difference between great-circle azimuths and
% rhumb azimuths. However, the differences in distances and azimuths are
% small (fractions of a degree in azimuth and little difference in the
% distance to the horizon) compared to the values from the loops
% (holdDistance and mean azimuth) in horz1fPts even when the DEM spans
% several degrees of longitude and the distance to the horizon is 100+ km.
% Therefore the threshold in function horz1fPts above is set to 100 km, but
% feel free to change.

hpu = unique(Hpt);
horzAng = ang;
horzDis = dis;
horzAzm = azm;
ii = (1:length(Hpt)).';
useGeodetic = Hpt~=ii & dis>threshold;
for jj=1:length(hpu)
    k = Hpt==hpu(jj) & useGeodetic;
    [horzAzm(k),horzAng(k),horzDis(k)] = geodetic2aer(lat(k),lon(k),Z(k),...
        lat(hpu(jj)),lon(hpu(jj)),Z(hpu(jj)),E);
end
tzero = horzAng<0;
if any(tzero)
    horzAng(tzero) = 0;
end
end