function fice=makeIceMask(S,RR,glims_res)
% make a fractional image of glacier ice
%input: S - output from shapread of GLIMS polygons w/ M and Z data removed
% RR- mapcellref
% glims_res - multiple of pixel size for rasterizing glims, e.g. 5
gmask=true(RR.RasterSize);
%create a higher res mask with 5% larger extent to deal with edge problems


pixelsize=[RR.CellExtentInWorldX/glims_res, ...
    RR.CellExtentInWorldY/glims_res];

%extend image 5 pixels to deal w/ edges when reprojected
Xlimit=[RR.XWorldLimits(1)-RR.CellExtentInWorldX*5 ...
        RR.XWorldLimits(2)+RR.CellExtentInWorldX*5];
Ylimit=[RR.YWorldLimits(1)-RR.CellExtentInWorldY*5 ...
        RR.YWorldLimits(2)+RR.CellExtentInWorldY*5];

[gmaskBig,RasterReferenceB]=rasterReprojection(gmask,RR,...
     'outProj',RR.ProjectedCRS,'method','nearest','PixelSize',...
     pixelsize,'XLimit',Xlimit,'YLimit',Ylimit);
gmaskBig(:)=false;

X=RasterReferenceB.XWorldLimits;
Y=RasterReferenceB.YWorldLimits;
bbox=[X(1) Y(1); X(2) Y(1);X(2) Y(2);X(1) Y(2)];
[lat,lon]=projinv(RR.ProjectedCRS,bbox(:,1),bbox(:,2));
XV=[min(lon); max(lon)];
YV=[min(lat); max(lat)];
for i=1:length(S)
    xx=S(i).X;
    yy=S(i).Y;
%     t=isnan(xx) | isnan(yy);
%     xx=xx(~t);
%     yy=yy(~t);
    in=inpolygon(xx,yy,XV,YV);
    if any(in(:))
        [x,y]=projfwd(RR.ProjectedCRS,yy,xx);
        [r,c]=map2pix(RasterReferenceB,x,y);
        %build mask for each NaN separated polygon
        t=find(isnan(xx) | isnan(yy));
        for j=1:length(t)
            if j==1
               ind=1:(t(j)-1); 
            else
               ind=(t(j-1)+1):(t(j)-1);
            end
            mask=poly2mask(c(ind),r(ind),...
                size(gmaskBig,1),size(gmaskBig,2));
            gmaskBig(mask)=true;
        end
    end
end

fice=rasterReprojection(single(gmaskBig),RasterReferenceB,'rasterref',RR);
fice(fice<0.01)=0;