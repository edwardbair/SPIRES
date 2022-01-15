function [SForward,SBackward] = internalHorizonPts(angToRotate,lat,lon,Z,E,UseParallel)
%rotation when not multiple of 90 degrees
assert(mod(angToRotate,90)~=0,...
    'rotation %f is a multiple of 90, use internal90 instead',angToRotate)

origSize = size(lat);

% keep track of the points' coordinates
[r,c] = ndgrid(1:size(Z,1),1:size(Z,2));
idx = sub2ind(size(Z),r,c);

%rotate the inputs, including a mask
%(accounting for possibility that some lat is zero or negative so set zeros
%to -Inf as a placeholder)
latmask = lat;
latmask(latmask==0) = -Inf;
lat = imrotate(lat,angToRotate);
lon = imrotate(lon,angToRotate);
Z = imrotate(Z,angToRotate);
idxR = imrotate(idx,angToRotate);
% set Z boundaries outside box to NaN, mask
latmask = imrotate(latmask,angToRotate);
mask = latmask==0;
Z(mask) = NaN;
lat(mask) = NaN;
lon(mask) = NaN;
idxR(idxR<=0) = NaN;

% output arrays
holdH = nan(size(lat));
holdA = nan(size(lat,2),2); % for mean of azimuths along profile
holdDis = nan(size(lat));
backH = nan(size(holdH));
backA = nan(size(holdA));
backDis = nan(size(holdDis));

if UseParallel
    parfor k=1:size(lat,2)
        [horzAng,horzDis,azm,npts] = horizonAlongProfile(lat(:,k),...
            lon(:,k),Z(:,k),E);
        holdH(:,k) = horzAng;
        holdA(k,:) = [atan2d(mean(sind(azm),'omitnan'),mean(cosd(azm),'omitnan'))...
            nnz(~isnan(npts))];
        holdDis(:,k) = horzDis;
        [horzAng,horzDis,azm,npts] = horizonAlongProfile(flipud(lat(:,k)),...
            flipud(lon(:,k)),flipud(Z(:,k)),E);
        backH(:,k) = horzAng;
        backA(k,:) = [atan2d(mean(sind(azm),'omitnan'),mean(cosd(azm),'omitnan'))...
            nnz(~isnan(npts))];
        backDis(:,k) = horzDis;
    end
else
    for k=1:size(lat,2)
        [horzAng,horzDis,azm,npts] = horizonAlongProfile(lat(:,k),...
            lon(:,k),Z(:,k),E);
        holdH(:,k) = horzAng;
        holdA(k,:) = [atan2d(median(sind(azm),'omitnan'),median(cosd(azm),'omitnan'))...
            nnz(~isnan(npts))];
        holdDis(:,k) = horzDis;
        [horzAng,horzDis,azm,npts] = horizonAlongProfile(flipud(lat(:,k)),...
            flipud(lon(:,k)),flipud(Z(:,k)),E);
        backH(:,k) = horzAng;
        backA(k,:) = [atan2d(median(sind(azm),'omitnan'),median(cosd(azm),'omitnan'))...
            nnz(~isnan(npts))];
        backDis(:,k) = horzDis;
    end
end
backH = flipud(backH);
backDis = flipud(backDis);
backA = flipud(backA);

SForward = cleanupRotated(origSize,angToRotate,idxR,holdA,holdH,holdDis);
SBackward = cleanupRotated(origSize,angToRotate,idxR,backA,backH,backDis);
end