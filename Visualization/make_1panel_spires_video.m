function make_1panel_spires_video(vidname,X,mdates,rrb)
%make a reprojected video of single variable only
%inputt
%vidname - filename, e.g. 2001h08v05.mp4
%X - cube of variable
%mdates - matlab dates
%rrb - rasterref w/ ProjectedCRS field
f=VideoWriter(vidname);

f.FrameRate=10;
f.Quality=90;
open(f);

f1=figure('Position',[100   1   1900  800],'Color','w');
    set(f1,'toolbar','none');
    set(gca,'NextPlot','replaceChildren');
      cm=colormap(parula);
cm(1,:)=[0.5 0.5 0.5];

[x,y]=worldGrid(rrb);
[lat,lon]=projinv(rrb.ProjectedCRS,x,y);
lat_l=ceil(lat(1,1)):-4:floor(lat(end,1));
lon_l=ceil(lon(1,1)):4:floor(lon(1,end));
    
[x,y]=projfwd(rrb.ProjectedCRS,mean(lat_l)*ones(size(lon_l)),lon_l);
[clon,~]=worldToIntrinsic(rrb,x,y);
    
[x,y]=projfwd(rrb.ProjectedCRS,lat_l,mean(lon_l)*ones(size(lat_l)));
[~,rlat]=worldToIntrinsic(rrb,x,y);
    imagesc;    
    colormap(cm);
    set(gca,'YDir','reverse','Box','on','XTicklabel',[],'YTickLabel',[]);
    set(gca,'YTick',rlat,'YTickLabel',num2str(lat_l'));
    set(gca,'XTick',clon,'XTickLabel',num2str(lon_l'))
for i=1:size(X,3)
    imagesc(X(:,:,i),'AlphaData',~isnan(X(:,:,i)));
    axis image;
    title(datestr(mdates(i)));
    frame=getframe(f1);
    writeVideo(f,frame)
end
close(f);
close(f1);
end