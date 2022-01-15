function [A,H,D] = horizonAllDirections(Z,R,varargin)
% [A,H,D] = horizonAllDirections(Z,R [,name-value pairs])
%calculates horizons for elevation grid in all azimuths around circle,
%either -180 to +180 degrees counter-clockwise or 0 to 360 degrees clockwise
%
%Values returned are the horizon angles and the azimuths for each point,
%along with the distances (meters) from each point to its horizon.
%
%Input
%   Z - elevation grid in meters, single or double precision
%   R - raster reference for the elevation grid, can be a MapCellsReference
%       (or MapPostingsReference) object, in which case the grid is projected
%       and distances between points are based on the projected coordinates,
%       or a GeographicCellsReference (or GeographicPostingsReference object),
%       in which case the distances between points are calculated from the
%       great-circle distances
%       If on MATLAB version R2020b or later, the R object should include
%       the ProjectedCRS or GeographicCRS field but if omitted or empty the
%       'proj' or 'planet' optional input can be used
%Optional input, name-value pairs, case-insensitive
%   'proj' - projection structure or projcrs object, needed if and only if
%       R is a MapCellsReference or MapPostingsReference object without the
%       non-empty ProjectedCRS field
%   'planet' - applicable only if R is a GeographicCellsReference or
%       GeographicPostingsReference object, default 'earth', value is
%       needed only if the R object does not contain a non-empty
%       GeographicCRS field
%   'nHorz' - number of horizon directions, default 64, but adjusted
%       to multiple of 4 otherwise (this number might be automatically
%       incremented by 1 or 2 to cover the entire azimuth circle)
%   'parallel' - followed by 'profile' to parallelize along the columns of
%       the rotated grid or by 'rotate' to run the suite of rotation angles
%       in parallel, default '' in which case no parallel processing is
%       employed
% Output
%   A - vector of azimuths, 1 for each rotation angle in nHorz directions
%   H - 3-D horizons in BSQ (band-sequential) format, the 3rd dimension
%       specifying the elevation angle to the horizon (i.e. elevation angle
%       above the horizontal)
%   D - distances to the horizon in meters corresponding to output H
%
% NOTE
% The 3D output is in 'bsq' (band-sequential) format. You can use bsq2bip
% and bip2bsq to switch between band-sequential and band-interleaved-by-line
% formats.
% The way azimuths are handled assumes that the origin of the grid is the
% NW corner. Depending on the return value from azimuthPreference, the
% azimuths either range over ±180° counter-clockwise with 0° south, or over
% 0° to 360° clockwise with 0° north.

% defaults
defaultH = 64;
defaultPlanet = 'wgs84';
minargs = 2;
maxargs = 10;
narginchk(minargs,maxargs)

p = inputParser;
addRequired(p,'Z',@(x) isnumeric(x) && ismatrix(x) && isfloat(x))
addRequired(p,'R',@(x) contains(class(x),'rasterref'))
addParameter(p,'nhorz',defaultH,@(x) isnumeric(x) && isscalar(x) && x>0)
addParameter(p,'planet',defaultPlanet,@ischar)
addParameter(p,'proj',struct([]),@(x) isstruct(x) ||...
    contains(class(x),'projcrs'))
addParameter(p,'parallel','',@ischar)
parse(p,Z,R,varargin{:});

% set parallel option
useParallelAlongProfiles = false;
useParallelAcrossRotations = false;
if strncmpi(p.Results.parallel,'pro',3)
    useParallelAlongProfiles = true;
elseif strncmpi(p.Results.parallel,'rot',3)
    useParallelAcrossRotations = true;
elseif ~isempty(p.Results.parallel)
    error('''parallel'' option ''%s'' not recognized',p.Results.parallel)
end

%other inputs
nHorz = p.Results.nhorz;
if contains(class(R),'geographic','IgnoreCase',true)
    if ~isempty(p.Results.proj)
        warning('R (2nd argument) is of class %s, so ''proj'' should be empty, it''s ignored',class(R))
    end
    proj = struct([]);
    useLatLon = true;
    if any(contains(fieldnames(R),'GeographicCRS'))
        if isempty(R.GeographicCRS) || verLessThan('map','5.0')
            E = referenceEllipsoid(p.Results.planet);
        else
            E = R.GeographicCRS.Spheroid;
        end
    else
        E = referenceEllipsoid(p.Results.planet);
    end
else
    assert(contains(class(R),'map.rasterref.Map'),...
        'R (2nd argument) must be a Map or Geographic raster reference')
    if any(contains(fieldnames(R),'ProjectedCRS'))
        if isempty(R.ProjectedCRS) || verLessThan('map','5.0')
            assert(~isempty(p.Results.proj),...
                'R (2nd argument) is of class %s, but ProjectedCRS is empty so ''proj'' can''t be empty')
            proj = p.Results.proj;
        else
            proj = R.ProjectedCRS;
        end
    else
        assert(~isempty(p.Results.proj),...
            'R (2nd argument) is of class %s, so ''proj'' can''t be empty',...
            class(R))
        proj = p.Results.proj;
    end
    useLatLon = false;
end

% nHorz must be multiple of 4
if nHorz<4
    nHorz = 4;
    warning('''nHorz'' must be >=4, set to 4')
elseif mod(nHorz,4)~=0
    nHorz = 4*ceil(nHorz/4);
    warning('''nHorz'' must be a multiple of 4, set to %d',nHorz)
end

% rotation angles, need just southern half of directions because horz2d
% computes forward and backward
ang = linspace(-90,90,1+nHorz/2);
ang(1) = []; % 90 and -90 gotten on same

% 3D arays to hold horizons in forward and backward directions for each azimuth
HgridF = zeros(size(Z,1),size(Z,2),length(ang),'single');
HgridB = zeros(size(HgridF),'like',HgridF);
AvecF = zeros(length(nHorz),1);
DgridF = zeros(size(HgridF),'like',HgridF);
DgridB = zeros(size(HgridF),'like',HgridF);
AvecB = zeros(length(nHorz),1);

% compute, depending on parallel option
if useParallelAcrossRotations
    % useParallelAlongProfiles passed to the rotational parfor loops, even
    % though this variable is set to false, but could be true in situations
    % where nested parallelism could be implemented
    if useLatLon
        passVals = {Z,R,useParallelAlongProfiles,E};
    else
        passVals = {Z,R,useParallelAlongProfiles,proj};
    end
    W = parallel.pool.Constant(passVals);
    if useLatLon
        parfor h=1:length(ang)
            pv = W.Value;
            [SF,SB] = horizonRotatedLatLon(ang(h),pv{1},pv{2},pv{3},pv{4});
            HgridF(:,:,h) = single(SF.horzAng);
            HgridB(:,:,h) = single(SB.horzAng);
            DgridF(:,:,h) = single(SF.horzDis);
            DgridB(:,:,h) = single(SB.horzDis);
            AvecF(h) = SF.azm;
            AvecB(h) = SB.azm;
        end
    else
        parfor h=1:length(ang)
            pv = W.Value;
            [SF,SB] = horizonRotatedProj(ang(h),pv{1},pv{2},pv{3},pv{4});
            HgridF(:,:,h) = single(SF.horzAng);
            HgridB(:,:,h) = single(SB.horzAng);
            DgridF(:,:,h) = single(SF.horzDis);
            DgridB(:,:,h) = single(SB.horzDis);
            AvecF(h) = SF.azm;
            AvecB(h) = SB.azm;
        end
    end
else
    if useLatLon
        for h=1:length(ang)
            [SF,SB] = horizonRotatedLatLon(ang(h),Z,R,useParallelAlongProfiles,E);
            HgridF(:,:,h) = single(SF.horzAng);
            HgridB(:,:,h) = single(SB.horzAng);
            DgridF(:,:,h) = single(SF.horzDis);
            DgridB(:,:,h) = single(SB.horzDis);
            AvecF(h) = SF.azm;
            AvecB(h) = SB.azm;
        end
    else
        for h=1:length(ang)
            [SF,SB] = horizonRotatedProj(ang(h),Z,R,useParallelAlongProfiles,proj);
            HgridF(:,:,h) = single(SF.horzAng);
            HgridB(:,:,h) = single(SB.horzAng);
            DgridF(:,:,h) = single(SF.horzDis);
            DgridB(:,:,h) = single(SB.horzDis);
            AvecF(h) = SF.azm;
            AvecB(h) = SB.azm;
        end
    end
end

%concatenate and sort by azimuths
Hcat = cat(3,HgridF,HgridB);
Dcat = cat(3,DgridF,DgridB);
Acat = cat(2,AvecF,AvecB);
[A,I] = sort(Acat);
H = zeros(size(Hcat),'like',Hcat);
D = zeros(size(Dcat),'like',Dcat);
for k=1:size(H,3)
    H(:,:,k) = Hcat(:,:,I(k));
    D(:,:,k) = Dcat(:,:,I(k));
end
mind = prctile(D(:),.5);
D(D<mind) = mind;

% smooth 3D H & D
kernelSize = round(max(3,length(A)/(128/5))); % 3x3x3 minimum, 5x5x5 if nHorz=128
H = smooth3(H,'gaussian',kernelSize);
D = smooth3(D,'gaussian',kernelSize);

% make periodic and return
[A,addFirst,addLast] = fillAzimuthEnds(A);
H = makeVariablePeriodic(addFirst,addLast,H,A);
D = makeVariablePeriodic(addFirst,addLast,D,A);
end

function [newA,addFirst,addLast] = fillAzimuthEnds(azm)
if min(azm)<0
    azrange = [-180 180];
else
    azrange = [0 360];
end
% add at beginning and/or end
addFirst = false;
addLast = false;
azm = azm(:).';
if min(azm)>min(azrange)
    azm = [min(azrange) azm];
    addFirst = true;
end
if max(azm)<max(azrange)
    azm = [azm max(azrange)];
    addLast = true;
end
newA = unique(azm); %should already, but just in case
end

function newX = makeVariablePeriodic(addFirst,addLast,X,azm)
% do nothing if already covers ends
if ~(addFirst || addLast)
    newX = X;
else
    X = bsq2bip(X);
    if length(azm)==size(X,1) % shouldn't, but just in case
        newX = bip2bsq(X);
    else
        newX = zeros(length(azm),size(X,2),size(X,3),'like',X);
        if addFirst && addLast
            newX(2:end-1,:,:) = X;
            newX(1,:,:) = (newX(2,:,:)+newX(end-1,:,:))/2;
            newX(end,:,:) = newX(1,:,:);
        elseif addFirst
            newX(2:end,:,:) = X;
            newX(1,:,:) = newX(end,:,:);
        else %addLast
            newX(1:end-1,:,:) = X;
            newX(end,:,:) = newX(1,:,:);
        end
        newX = bip2bsq(newX);
    end
end
end