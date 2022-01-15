function vf = viewFactor(A,H,Z,R,varargin)
% vf = viewFactor(A,H,Z,R)
% calculate view factor for topographic grid after running horizonAllDirections
%
% Input
%   A,H - output from horizonAllDirections
%       A - vector of azimuths at which horizons are calculated
%       H - 3-D horizons in BSQ (band-sequential) format, the 3rd dimension
%       specifying the elevation angle to the horizon (i.e. elevation angle
%       above the horizontal)
%   Z - 2-D elevation grid in meters, single or double precision
%       (dimensions must match the first two dimensions of the H cube)
%   R - raster reference for the elevation grid, can be a MapCellsReference
%       (or MapPostingsReference) object, in which case the grid is projected
%       and distances between points are based on the projected coordinates,
%       or a GeographicCellsReference (or GeographicPostingsReference object),
%       in which case the distances between points are calculated from the
%       great-circle distances
% Optional input (just a scalar, not a name-value pair)
%   useParallel - true (default) or false (or 1 or 0)
%
% Output
%   vf - sky view factor for the input grid

p = inputParser;
addRequired(p,'A',@(x) isvector(x) && isnumeric(x))
addRequired(p,'H',@(x) ndims(x)==3 && isnumeric(x))
addRequired(p,'Z',@(x) isnumeric(x) && ismatrix(x) && isfloat(x))
addRequired(p,'R',@(x) contains(class(x),'rasterref'))
addOptional(p,'useParallel',true,@(x) isscalar(x) && (islogical(x) || isnumeric(x)))
parse(p,A,H,Z,R,varargin{:})

useParallel = logical(p.Results.useParallel);

% check dimensions
assert(isequal(size(Z),[size(H,1) size(H,2)]),...
    'size of elevation grid Z [%d %d] must be same as first two dimensions of horizon cube H [%d %d]',...
    size(Z,1),size(Z,2),size(H,1),size(H,2))
assert(length(A)==size(H,3),'length of A vector [%d] must be same as 3rd dimension of horizon cube H [%d]',...
    length(A),size(H,3))

% slope and aspect of topographic grid
[slope,aspect] = topographicSlope(Z,R);

% convert band-sequential H to band-interleaved-by-pixel so that horizon
% angles for a point are together
H = bsq2bip(H);

%flatten slope, aspect, and horizons to run the loop in just 1 dimension
slope = slope(:);
aspect = aspect(:);
H = reshape(H,size(H,1),size(H,2)*size(H,3));

%output
vf = zeros(size(slope));

% loop through all points
if useParallel
    parfor k=1:length(slope)
        vf(k) = viewFpt(A,H(:,k),slope(k),aspect(k));
    end
else
    for k=1:length(slope)
        vf(k) = viewFpt(A,H(:,k),slope(k),aspect(k));
    end
end

% convert output to matrix
vf = reshape(vf,size(Z));

end

function V = viewFpt(azm,horzAng,slopeDegrees,aspectDegrees )
% [ V ] = viewFpt(azm,horzAng,slopeDegrees,aspectDegrees )
%
%sky view factor for horizon circle for a point
% (fraction of the sky open to that point in grid)
%based on equation 7b in Dozier, J., and J. Frew (1990), Rapid calculation
%   of terrain parameters for radiation modeling from digital elevation data,
%   IEEE Trans. Geosci. Remote Sens., 28, 963-969, doi: 10.1109/36.58986.
%
%INPUT
%   azm - azimuths of horizon vector in degrees
%   horzAng - horizon vector (elevation angle from horizontal) in degrees
%   slopeDegrees - slope from horizontal
%   aspectDegrees - Angles for azm and aspect can be either 0° to 360°
%       (clockwise from north) or ±180° (counter-clockwise, 0° south).
%       Either can be used but make sure that aspect and azm use the same
%       coordinate system.
%
%OUTPUT
%   V - view factor for this point

% expand azimuths to cover full circle, whether ±180° or 0° to 360° and
% convert to row vectors
azm = azm(:).';
horzAng = horzAng(:).';

% range check
if min(azm)<0 % in ±180° range
    assert(min(azm)>=-180 && max(azm)<=180,'azimuths must lie within +/-180')
    if min(azm)>-180
        azm = [-180 azm];
        horzAng = [horzAng(1) horzAng];
    end
    if max(azm)<180
        azm = [azm 180];
        horzAng = [horzAng horzAng(end)];
    end
else % in 0° to 360° range
    assert(min(azm)>=0 && max(azm)<=360,'azimuths must lie within 0 to 360')
    if min(azm)>0
        azm = [0 azm];
        horzAng = [horzAng(1) horzAng];
    end
    if max(azm)<360
        azm = [azm 360];
        horzAng = [horzAng horzAng(end)];
    end
end

azmRadian = (pi/180)*azm;
horzAng(horzAng<0) = 0;

% convert output from horizon program, radians and from zenith
% integrand equivalent eq 7b in Dozier & Frew 1990
H = (pi/180)*(90-horzAng(:).');
if slopeDegrees>0
    aspectRadian = (pi/180)*(aspectDegrees);
    qIntegrand = cosd(slopeDegrees)*(1-cos(2*H))+...
        cos(azmRadian-aspectRadian)*sind(slopeDegrees).*(2*H-sin(2*H));
else
    qIntegrand = 1-cos(2*H);
end
V = trapz(azmRadian,qIntegrand)/(4*pi);
end