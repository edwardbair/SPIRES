function [horzangles,azm ] = horz1dForward( lat,lon,Z,varargin )
% [horzangles,azm ] = horz1dForward( lat,lon,Z, {planet} )
%calculates horizon points and angles in the forward direction along an
%elevation profile
% (to go in the backward direction, just flip the input vectors)
%
% based on - Dozier, J., J. Bruno, and P. Downey (1981), A faster solution
% to the horizon problem, Comp. Geosci., 7, 145-151,
% doi: 10.1016/0098-3004(81)90026-1
% adapted here to lat-lon coordinates and accounting for planet curvature
%
% INPUT - vectors of same length, need not be equispaced
%   lat - latitude vector of points
%   lon - longitude vector
%   Z - elevations
% OPTIONAL INPUT
%   planet - name of planet as a character string, case insensitive
%
% OUTPUT
%   horzangles - angles to the horizon in the forward direction
%   azm - azimuths to the horizon in the forward direction (i.e.
%       corresponding to the sun's azimuth)

% check input
assert(isvector(lat) && isvector(lon) && isvector(Z),...
    'input lat, lon, and Z must be vectors')
assert(isequal(size(lat),size(lon),size(Z)),...
    'input lat, lon, and Z must be vectors of same size')
if ~isfloat(Z)
    Z = double(Z);
end

optarg = nargin-3;
if optarg ==1
    S = referenceSphere(varargin{1});
else
    S = referenceSphere('earth');
end
R = S.Radius;

phi = distance(lat,lon,lat(1),lon(1)); % arc distances - assume spherical
[horzpts, horzangles] = horz1f(phi,Z,R);
% azimuths to each point's horizon
% (nb: if we ever want to save distances to the horizons, we can do it
% here)
azm = nan(size(horzpts));
for k=1:length(azm)
    if horzpts(k)~=k % leave azm as NaN when point is its own horizon
        azm(k) = azimuth(lat(k),lon(k),lat(horzpts(k)),lon(horzpts(k)));
    end
end
azm = 180-azm;
end

function [H,Ha] = horz1f(phi,Z,R)

H = zeros(size(phi));
Ha = zeros(size(phi));

% end point is its own horizon
N = length(H);
H(end) = N;

for i=N-1:-1:1 % loop from next-to-backward to beginning
    zi = Z(i);
    phii = phi(i);
    % Start with next-to-adjacent point. We don't consider the adjacent
    % point itself; this seems to reduce noise.
    k = i+2;
    if k > N
        k = N;
    end
    
    % loop until horizon found - note the use of flag to simulate the
    % 'do while' construct in C
    flag = true;
    while flag
        j = k;
        k = H(j);
        sij = hangle(phii,phi(j),zi,Z(j),R);
        sihj = hangle(phii,phi(k),zi,Z(k),R);
        % if sij>=sihj, horizon has been found, otherwise set j to its
        % horizon and loop again
        if sij>=sihj
            flag = false;
        end
    end
    % if angle(i,j)>angle(i,H(j)), j is i's horizon, but if angle
    % is <=0, i is its own horizon, otherwise angle(i,j)=angle(i,H(j))
    if sij>sihj
        H(i) = j;
        Ha(i) = sij;
    elseif sij==0
        H(i) = i;
        Ha(i) = sij;
    else
        H(i) = k;
        Ha(i) = sihj;
    end
end
end

function H = hangle(phi1,phi2,z1,z2,R)
if z2<=z1
    H = 0;
else
    Z = z2-z1;
    phi = abs(phi2-phi1);
    % horizon angle
    H = max(0,rad2deg(atan2(cosd(phi)*(R+Z)-R,sind(phi)*R+Z)));
end
end