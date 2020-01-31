function [ hForward, hBackward, azmF, azmB ] = horz2d( lat,lon,Z,azm,varargin )
% [ hForward, hBackward, azmF, azmB ] = horz2d( lat,lon,Z,azm, {planet} )
%forward and backward horizon angles over elevation grid, where azimuths
%can vary along columns (especially with projected elevations)
%
%INPUT
%   lat,lon - grids of geographic coordinates in degrees
%   Z - elevation grid in meters
%   azm - scalar, in degrees, 0 south, + counter-clockwise, measured
%       with respect to the grid itself, which may not be wrt north
%
% OPTIONAL INPUT
%   planet - name of planet as a character string, case insensitive
%
%OUTPUT
% horizons that contains horizon angles for 2 directions, forward & backward
%   hForward - horizons in forward direction (direction azmF)
%   hBackward - horizons in backward direction (direction azmB)
% Azimuths - azimuth directions, size(Z), 0 south, + ccw (the Backward
%   direction is +/- 180 deg from the forward
%
% NOTE
% The way azimuths are handled assumes that the origin of the grid is the
% NW corner

% check inputs
assert(isequal(size(lat),size(lon),size(Z)),...
    'input lat, lon, and Z must be 2D and of equal size')
assert(isscalar(azm),...
    'input azimuth must be a scalar')

optarg = nargin-4;
if optarg==1
    planet = varargin{1};
else
    planet = 'earth';
end

% method is to rotate the elevation, lat, & lon grids, compute the
% horizons column-wise, and then rotate the results back

% rotate the input grids, with rback & mask needed to rotate back later
[ZR,rback,mask] = RotateHorz(Z,azm);
lat = RotateHorz(lat,azm,'nomask');
lon = RotateHorz(lon,azm,'nomask');

Hf = zeros(size(ZR));
Hb = zeros(size(Hf));
Af = zeros(size(Hf));
Ab = zeros(size(Hb));

% special case, multiple of 90 deg
if rback.multipleof90
    parfor c=1:size(ZR,2) % forward direction
        [Hf(:,c),Af(:,c)] = horz1dForward(lat(:,c),lon(:,c),ZR(:,c),planet);
    end
    % flip for the backward direction
    lat = flipud(lat);
    lon = flipud(lon);
    ZR = flipud(ZR);
    parfor c=1:size(ZR,2) % backward direction
        [Hb(:,c),Ab(:,c)] = horz1dForward(lat(:,c),lon(:,c),ZR(:,c),planet);
    end
    Hb = flipud(Hb);
    Ab = flipud(Ab);
    % rotate the results back
    rback.minmax = [min(Hf(:)) max(Hf(:))];
    horzAng = RotateHorz(Hf,rback);
    rback.minmax = [min(Hb(:)) max(Hb(:))];
    horzB = RotateHorz(Hb,rback);
    rback.minmax = [min(Af(:)) max(Af(:))];
    azm1 = RotateHorz(Af,rback);
    rback.minmax = [min(Ab(:)) max(Ab(:))];
    azm2 = RotateHorz(Ab,rback);
else
    parfor c=1:size(ZR,2)
        clat = lat(:,c);
        clon = lon(:,c);
        cZ = ZR(:,c);
        t = mask(:,c);
        cH = Hf(:,c);
        cA = Af(:,c);
        if nnz(t)>1
            [cH(t),cA(t)] = horz1dForward(clat(t),clon(t),cZ(t),planet);
        end
        Hf(:,c) = cH;
        Af(:,c) = cA;
    end
    % flip for the backward direction
    lat = flipud(lat);
    lon = flipud(lon);
    ZR = flipud(ZR);
    mask = flipud(mask);
    parfor c=1:size(ZR,2)
        clat = lat(:,c);
        clon = lon(:,c);
        cZ = ZR(:,c);
        t = mask(:,c);
        cH = Hb(:,c);
        cA = Ab(:,c);
        if nnz(t)
            [cH(t),cA(t)] = horz1dForward(clat(t),clon(t),cZ(t),planet);
        end
        Hb(:,c) = cH;
        Ab(:,c) = cA;
    end
    Hb = flipud(Hb);
    Ab = flipud(Ab);
    % rotate back, but minmax needs resetting
    rback.minmax = [min(Hf(:)) max(Hf(:))];
    horzAng = RotateHorz(Hf,rback);
    rback.minmax = [min(Hb(:)) max(Hb(:))];
    horzB = RotateHorz(Hb,rback);
    rback.minmax = [min(Af(:)) max(Af(:))];
    azm1 = RotateHorz(Af,rback);
    rback.minmax = [min(Ab(:)) max(Ab(:))];
    azm2 = RotateHorz(Ab,rback);
end

% eliminate bad values and set azimuths for points who are their own
% horizons to neighboring azimuths

mask1 = horzAng<0 | isnan(horzAng) | isnan(azm1);
mask2 = horzB<0 | isnan(horzB) | isnan(azm2);
azm1(mask1) = 0;
azm2(mask2) = 0;
azm1 = inpaintCoherent(azm1,mask1);
azm2 = inpaintCoherent(azm2,mask2);

azmF = azm1;
azmB = azm2;
hForward = horzAng;
hBackward = horzB;

end