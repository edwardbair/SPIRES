function [ newInR, newRaster ] = coarsenInput(raster,InRR,OutRR,planet,method)
%coarsen input raster if resolution of output is significantly coarser
%Input
%   raster - original raster
%   InRR - input raster reference
%   OutRR - output raster reference
%   planet - designation of which planet
%   method - interpolation method
%Output
%   newInR - raster reference of newRaster
%   newRaster - coarsened input raster

thresholdRatio = 1.2;
[x11,y11,dx,dy,inProj] = cornerCoords(InRR);
[~,~,odx,ody,outProj] = cornerCoords(OutRR);
hdy = dy;
hdx = dx;
hody = ody;
hodx = odx;

if xor(inProj,outProj) % one projected, one geographic
    S = referenceEllipsoid(planet);
    R = S.MeanRadius;
    % convert the geographic one to projected
    if inProj
        hody = deg2rad(ody)*R;
        hodx = deg2rad(odx)*R*(cosd(mean(OutRR.LatitudeLimits)));
    else
        hdy = deg2rad(dy)*R;
        hdx = deg2rad(dx)*R*(cosd(mean(InRR.LatitudeLimits)));
    end
end
minRatio = min(abs([hody hodx]./[hdy hdx]));
iterations = floor(log2(minRatio/thresholdRatio));

% reduce by about half every interation
for k=1:iterations
    if strcmpi(method,'nearest')
        raster = medfilt2(raster,'symmetric');
        newRaster = imresize(raster,1/2,method);
    else
        newRaster = impyramid(raster,'reduce');
        % use nearest neighbor to keep NaNs from propagating
        if any(isnan(newRaster(:)))
            X = imresize(raster,[size(newRaster,1) size(newRaster,2)],'nearest');
            t = isnan(newRaster);
            newRaster(t) = X(t);
        end
    end
    adjustment = size(raster)./size(newRaster);
    newdy = dy*adjustment(1);
    newdx = dx*adjustment(2);
    x11 = x11+(newdx-dx);
    y11 = y11+(newdy-dy);
    raster = newRaster;
    dx = newdx;
    dy = newdy;
end

% adjust raster reference if we did anything
if iterations>=1
    RM = makerefmat(x11,y11,dx,dy);
    if inProj
        newInR = refmatToMapRasterReference(RM,size(newRaster));
    else
        newInR = refmatToGeoRasterReference(RM,size(newRaster));
    end
else
    newInR = InRR;
    newRaster = raster;
end
end