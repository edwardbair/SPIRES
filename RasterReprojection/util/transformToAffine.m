function [newRR] = transformToAffine(RR,rotationAngle)
%transform raster reference to same with rotation, an affine transform
cx = [1 2 2 1];
cy = [1 1 2 2];
if contains(class(RR),'geographic','IgnoreCase',true)
    [lat,lon] = intrinsicToGeographic(RR,cx,cy);
    quadin = polyshape(lon,lat);
    quadout = rotate(quadin,rotationAngle,[mean(lon) mean(lat)]);
    dydy = quadout.Vertices(4,2)-quadout.Vertices(1,2);
    dxdy = quadout.Vertices(4,1)-quadout.Vertices(1,1);
    dydx = quadout.Vertices(2,2)-quadout.Vertices(1,2);
    dxdx = quadout.Vertices(2,1)-quadout.Vertices(1,1);
    dx = [dxdy dxdx];
    dy = [dydy dydx];
    RM = makerefmat(lon(1),lat(1),dx,dy);
    newRR = refmatToGeoRasterReference(RM,RR.RasterSize,...
        RR.RasterInterpretation);
elseif contains(class(RR),'mapcells','IgnoreCase',true) ||...
        contains(class(RR),'mappostings','IgnoreCase',true)
    [x,y] = intrinsicToWorld(RR,cx,cy);
    quadin = polyshape(x,y);
    quadout = rotate(quadin,rotationAngle,[mean(x) mean(y)]);
    dydy = quadout.Vertices(4,2)-quadout.Vertices(1,2);
    dxdy = quadout.Vertices(4,1)-quadout.Vertices(1,1);
    dydx = quadout.Vertices(2,2)-quadout.Vertices(1,2);
    dxdx = quadout.Vertices(2,1)-quadout.Vertices(1,1);
    dx = [dxdy dxdx];
    dy = [dydy dydx];
    RM = makerefmat(x(1),y(1),dx,dy);
    newRR = refmatToMapRasterReference(RM,RR.RasterSize,...
        RR.RasterInterpretation);
else
    error('class input raster reference %s not recognized', class(RR))
end