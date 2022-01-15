function [tmpS] = cleanupRotated(origSize,rotatedAngle,idxR,A,H,D)
% clean up, remove NaNs because some coordinates are missed in the
% rotation, and some missed horizon angles are zero
origTrimMask = isnan(H);
% find the coordinates of the rotated and re-rotated image corners
gridID = true(origSize);
gridR = imrotate(gridID,rotatedAngle);
assert(isequal(size(gridR),size(H)),'rotated grid sizes not the same')
gridRback = imrotate(gridR,-rotatedAngle,'crop');
% find coordinates of corners region within gridRback (approximately a
%   quadrilateral)
[r,c] = find(gridRback);
rf = min(r);
rl = max(r);
cf = min(c);
cl = max(c);
for p=1:2
    Destination = nan(origSize);
    switch p
        case 1
            Source = H;
            origin = ~origTrimMask;
            minValue = nanmin(H(H>0));
            reRotated = imrotate(H,-rotatedAngle,'crop');
        case 2
            % distances at the end of the grid could be zero
            Source = D;
            origin = ~origTrimMask;
            minValue = nanmin(D(D>0));
            reRotated = imrotate(D,-rotatedAngle,'crop');
    end
    % size of re-rotated grid should be same
    reRotated = reRotated(rf:rl,cf:cl);
    if isempty(reRotated)
        pause
    end
    if ~isequal(size(reRotated),origSize)
        reRotated = imresize(reRotated,origSize,'nearest');
    end
    % fill the rotated values by indexed coordinates, then fill
    % their NaNs with the re-rotated values that are not NaN
    dest = idxR(origin);
    Destination(dest) = Source(origin);
    t = isnan(Destination);
    Destination(t) = reRotated(t);
    % NaNs exist because some coordinates are missed in the rotation, but
    % the matrices were initialized with NaNs so can fill
    t = isnan(Destination);
    if any(t,'all')
        X = Destination;
        X(t) = -1;
        % preserve zeros
        tsame = X==0;
        % fill in NaNs
        Destination = inpaintCoherent(X,t);
        % some at edges or within not filled by inpaintCoherent, use regionfill
        % for those
        t = isnan(Destination) | (Destination<minValue & ~tsame);
        if any(t,'all')
            Destination = regionfill(Destination,t);
        end
    end
    switch p
        case 1
            tmpS.horzAng = Destination;
        case 2
            tmpS.horzDis = Destination;
    end
    
    % weighted mean of the columwise azimuths
    wt = A(:,2)/nansum(A(:,2));
    azmMean = atan2d(nansum(sind(A(:,1)).*wt),...
        nansum(cosd(A(:,1)).*wt));
    % atan2d returns +/- 180, fix if we use 0 to 360
    if ~azimuthPreference && azmMean<0
        azmMean = 360+azmMean;
    end
    tmpS.azm = azmMean;
end
end