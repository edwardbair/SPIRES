function [ newRaster, newInRR ] = coarsenInput(raster,InRR,OutRR,planet,fillvalue)
%coarsen input raster if resolution of output is significantly coarser
%Input
%   raster - original raster
%   InRR - input raster reference
%   OutRR - output raster reference
%   planet - designation of which planet
%Output
%   newInR - raster reference of newRaster
%   newRaster - coarsened input raster

thresholdRatio = 3;
inProj = ~contains(class(InRR),'geographic','IgnoreCase',true);
outProj = ~contains(class(OutRR),'geographic','IgnoreCase',true);

if xor(inProj,outProj) % one projected, one geographic
    S = referenceEllipsoid(planet);
    R = S.MeanRadius;
    % convert the geographic pixelsize to meters
    if inProj
        hody = deg2rad(OutRR.CellExtentInLatitude)*R;
        hodx = deg2rad(OutRR.CellExtentInLongitude)*R*(cosd(mean(OutRR.LatitudeLimits)));
        hdy = InRR.CellExtentInWorldY;
        hdx = InRR.CellExtentInWorldX;
    else
        hdy = deg2rad(InRR.CellExtentInLatitude)*R;
        hdx = deg2rad(InRR.CellExtentInLongitude)*R*(cosd(mean(InRR.LatitudeLimits)));
        hody = OutRR.CellExtentInWorldY;
        hodx = OutRR.CellExtentInWorldX;
    end
else
    % either both geographic or both projected
    if inProj
        hody = OutRR.CellExtentInWorldY;
        hodx = OutRR.CellExtentInWorldX;
        hdy = InRR.CellExtentInWorldY;
        hdx = InRR.CellExtentInWorldX;
    else
        hody = OutRR.CellExtentInLatitude;
        hodx = OutRR.CellExtentInLongitude;
        hdy = InRR.CellExtentInLatitude;
        hdx = InRR.CellExtentInLongitude;
    end
end

cellRatio = mean([hodx/hdx hody/hdy]);

%use mapresize or georesize, depending on whether input is projected or
%geographic
if cellRatio>thresholdRatio
    if iscategorical(raster)
        method = 'nearest';
    else
        method = 'cubic';
    end
    A = double(raster);
    if ~isfloat(raster)
        A(raster==fillvalue) = NaN;
    end
    if inProj
        [newRaster,newInRR] = mapresize(A,InRR,2/cellRatio,method);
        if ~strcmpi(method,'nearest') && any(isnan(newRaster),'all')
            % use nearest neighbor to keep NaNs from propagating
            t = isnan(newRaster);
            [holdRaster,~] = mapresize(A,InRR,2/cellRatio,'nearest');
            newRaster(t) = holdRaster(t);
        end
    else
        [newRaster,newInRR] = georesize(A,InRR,2/cellRatio);
        if ~strcmpi(method,'nearest') && any(isnan(newRaster),'all')
            % use nearest neighbor to keep NaNs from propagating
            t = isnan(newRaster);
            [holdRaster,~] = georesize(A,InRR,2/cellRatio,'nearest');
            newRaster(t) = holdRaster(t);
        end
    end
    t = isnan(newRaster);
    if isa(raster,'single')
        newRaster = single(newRaster);
        if ~isnan(fillvalue)
            newRaster(t) = fillvalue;
        end
    elseif ~isa(raster,'double') && ~islogical(raster)
        newRaster = cast(round(newRaster),class(raster));
        newRaster(t) = fillvalue;
    elseif logical(raster)
        newRaster = newRaster>0.5;
        newRaster(t) = fillvalue;
    elseif ~isnan(fillvalue) % double is the only class left
        newRaster(t) = fillvalue;
    end
else
    newInRR = InRR;
    newRaster = raster;
end
end