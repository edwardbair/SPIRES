function [SForward,SBackward] = internal90(angToRotate,lat,lon,Z,E,UseParallel)
%straightforward rotation when multiple of 90 degrees
assert(mod(angToRotate,90)==0,...
    'this routine only applies to rotation angles that are a multiple of 90 degrees, use internalHorizonPts instead')
rk = round(angToRotate/90);

% if we want to keep track of the points
% [r,c] = ndgrid(1:size(Z,1),1:size(Z,2));
% idx = sub2ind(size(Z),r,c);

if rk~=0
    lat = rot90(lat,rk);
    lon = rot90(lon,rk);
    Z = rot90(Z,rk);
end

for iteration=1:2
    if iteration==2
        lat = flipud(lat);
        lon = flipud(lon);
        Z = flipud(Z);
    end
    % output arrays
    holdH = nan(size(lat));
    holdA = nan(size(lat,2),2); % for mean of azimuths along profile
    holdDis = nan(size(lat));
    
    if UseParallel
        parfor k=1:size(lat,2)
            [horzAng,horzDis,azm,horzPts] = horizonAlongProfile(lat(:,k),...
                lon(:,k),Z(:,k),E);
            holdH(:,k) = horzAng;
            holdA(k,:) = [atan2d(mean(sind(azm),'omitnan'),mean(cosd(azm),'omitnan'))...
                nnz(~isnan(horzPts))];
            holdDis(:,k) = horzDis;
        end
    else
        for k=1:size(lat,2)
            [horzAng,horzDis,azm,horzPts] = horizonAlongProfile(lat(:,k),...
                lon(:,k),Z(:,k),E);
            holdH(:,k) = horzAng;
            holdA(k,:) = [atan2d(mean(sind(azm),'omitnan'),mean(cosd(azm),'omitnan'))...
                nnz(~isnan(horzPts))];
            holdDis(:,k) = horzDis;
        end
    end
    
    if iteration==2
        holdH = flipud(holdH);
        holdDis = flipud(holdDis);
    end
    
    %rotate back
    if rk==0
        tmpS.horzAng = holdH;
        tmpS.horzDis = holdDis;
    else
        tmpS.horzAng = rot90(holdH,-rk);
        tmpS.horzDis = rot90(holdDis,-rk);
    end
    
    % weighted mean of the columwise azms
    wt = holdA(:,2)/sum(holdA(:,2),'omitnan');
    tmpS.azm = atan2d(sum(sind(holdA(:,1)).*wt,'omitnan'),...
        sum(cosd(holdA(:,1)).*wt,'omitnan'));
    % atan2d returns +/- 180, fix if we use 0 to 360
    if ~azimuthPreference && tmpS.azm<0
        tmpS.azm = 360+tmpS.azm;
    end
    
    % take care of any NaNs
    t = isnan(tmpS.horzAng);
    if any(t,'all')
        tmpS.horzAng(t) = 0;
        tmpS.horzAng = inpaintCoherent(tmpS.horzAng,t);
        % any NaN at edges
        t = isnan(tmpS.horzAng);
        if any(t,'all')
            tmpS.horzAng = regionfill(tmpS.horzAng,t);
        end
    end
    t = isnan(tmpS.horzDis);
    if any(t,'all')
        tmpS.horzDis(t) = 0;
        tmpS.horzDis = inpaintCoherent(tmpS.horzDis,t);
        % any NaN at edges
        t = isnan(tmpS.horzDis);
        if any(t,'all')
            tmpS.horzDis = regionfill(tmpS.horzDis,t);
        end
    end
    
    if iteration==2
        SBackward = tmpS;
    else
        SForward = tmpS;
    end
end
end