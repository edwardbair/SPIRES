function gfrac=makeIceMask(S,hdr,glims_res)
% make a fractional image of glacier ice
%input: S - output from shapread of GLIMS polygons w/ M and Z data removed
% hdr - hdr geo info for target
% glims_res - spatial resolution for rasterizing glims data, e.g. 30 m
gmask=false(hdr.RasterReference.RasterSize);
[gmaskBig,RefMatrixB,RasterReferenceB]=rasterReprojection(gmask,hdr.RefMatrix,hdr.ProjectionStructure,...
    hdr.ProjectionStructure,'PixelSize',[glims_res glims_res],'Method',...
    'nearest');

X=RasterReferenceB.XWorldLimits;
Y=RasterReferenceB.YWorldLimits;
bbox=[X(1) Y(1); X(2) Y(1);X(2) Y(2);X(1) Y(2)];
[lat,lon]=minvtran(hdr.ProjectionStructure,...
    bbox(:,1),bbox(:,2));
XV=[min(lon); max(lon)];
YV=[min(lat); max(lat)];
for i=1:length(S)
    xx=S(i).X;
    yy=S(i).Y;
    t=isnan(xx) | isnan(yy);
    xx=xx(~t);
    yy=yy(~t);
    in=inpolygon(xx,yy,XV,YV);
    if any(in(:))
        [x,y]=mfwdtran(hdr.ProjectionStructure,yy,xx);
        [r,c]=map2pix(RefMatrixB,x,y);
        mask=poly2mask(c,r,size(gmaskBig,1),size(gmaskBig,2));
        gmaskBig=gmaskBig | mask;
    end
end

gfrac=rasterReprojection(single(gmaskBig),RefMatrixB,...
    hdr.ProjectionStructure,hdr.ProjectionStructure,...
    'rasterref',hdr.RasterReference);
gfrac(gfrac<0.01)=0;