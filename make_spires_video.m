function make_spires_video(infiles,target,pshape)
%create reprojected spires MODIS video
%infiles - cell, N*1 list of h5 files to read
%output struct from smooth_and_run_modis
%w fields:
%fsca (canopy adj),
%grain radius (um)
%dust (ppmw)
%matdates (datenums)
%and hdr struct w it's own fields RefMatrix and ProjectionStructure
% target - hdr struct w. Projection Structure RasterReference to
%reproject to

%pshape - polyshape of boundary area
%topofile - h5 topofile location
fname='spires_video.avi';
f=VideoWriter(fname);

f.FrameRate=10;
f.Quality=90;
open(f);

dvisflag=false;
vars={'snow_fraction','grain_size','dust','albedo'};
info=h5info(infiles{1});
for i=1:length(info.Groups.Groups.Datasets)
    if strcmp(info.Groups.Groups.Datasets(i).Name,'deltavis')
        dvisflag=true;
        vars={'snow_fraction','grain_size','deltavis','albedo'};        
    end
end

f1=figure('Position',[100   1   700   900]);
set(f1,'toolbar','none');
ha=tight_subplot(2, 2, [0.025 0.001], [0.025 0.01], [0.025 0]);
set(ha,'NextPlot','replaceChildren');

[lon,lat]=pshape.boundary;

t=~isnan(lat) & ~isnan(lon);
[x,y]=mfwdtran(target.ProjectionStructure,lat,lon);
[row,col]=map2pix(target.RasterReference,x,y);

mask=poly2mask(col(t),row(t),target.RasterReference.RasterSize(1),...
    target.RasterReference.RasterSize(2));

[lon,lat]=pshape.boundingbox;
[x,y]=mfwdtran(target.ProjectionStructure,lat,lon);
[bbox_y,bbox_x]=map2pix(target.RasterReference,x,y);

cm=colormap(parula);

cm(1,:)=[0.5 0.5 0.5];

% [Slope,hdr]=GetTopography(topofile,'Slope');
% Aspect=GetTopography(topofile,'elevation');
% 
% Slope=rasterReprojection(Slope,hdr.RefMatrix,hdr.ProjectionStructure,....
%     target.ProjectionStructure,'rasterref',target.RasterReference);
% Aspect=rasterReprojection(Aspect,hdr.RefMatrix,hdr.ProjectionStructure,....
%     target.ProjectionStructure,'rasterref',target.RasterReference);

[x,y]=pixcenters(target.RefMatrix,target.RasterReference.RasterSize,'makegrid');
[lat,lon]=minvtran(target.ProjectionStructure,x,y);

for j=1:length(vars)
    axes(ha(j));
    ax=gca;
    imagesc;
    axis image;
    
    colormap(cm);
%     xlim([bbox_x(1) bbox_x(2)+70]);
%     ylim([bbox_y(2) bbox_y(1)+130]);
    
    
%      ax.XAxis.Color = 'w';
%      ax.YAxis.Color = 'w';
%     
    
    lat_l=ceil(lat(1,1)):-1:floor(lat(end,1));
    lon_l=ceil(lon(1,1)):1:floor(lon(1,end));
    
    [x,y]=mfwdtran(target.ProjectionStructure ,40*ones(size(lon_l)),lon_l);
    [~,clon]=map2pix(target.RefMatrix,x,y);
    
    [x,y]=mfwdtran(target.ProjectionStructure,lat_l,-120*ones(size(lat_l)));
    [rlat,~]=map2pix(target.RefMatrix,x,y);
    
   
    
    set(ax,'YDir','reverse','Box','on','XTick',clon,'XTicklabel',[],...
            'YTick',rlat,'YTickLabel',[]);
    
    if j==1 || j==3
        set(ax,'YTick',rlat,'YTickLabel',num2str(lat_l'));
    end
    if j>2
        set(ax,'XTick',clon,'XTickLabel',num2str(lon_l'))
    end
        
     ltr={'(a)','(b)','(c)','(d)'};
     if j==1
        c1=colorbar('Location','south');
        c1.Label.String='fsca';
%         c1.Label.Color=[1 1 1];
        caxis([0 1]);
     elseif j==2
        c2=colorbar('Location','south');
        c2.Label.String='grain radius, \mum';
%         c2.Label.Color=[1 1 1];
        caxis([40 900])
    elseif j==3
        c3=colorbar('Location','south');
        if dvisflag
            c3.Label.String='deltavis';
            caxis([0 0.4]);
        else
            c3.Label.String='dust conc, ppmw';
            caxis([0 1000]);
        end
%         c3.Label.Color=[1 1 1];
     elseif j==4
        c4=colorbar('Location','south');
        c4.Label.String='albedo';
%         c4.Label.Color=[1 1 1];
        caxis([0.4 0.9]);
     end
    
end

for ii=1:size(infiles,1)
    fname=infiles{ii};
    in=struct();
    for j=1:3
        if j==1
            [in.(vars{j}),in.matdates,in.hdr]=...
                GetEndmember(fname,vars{j});
        else
            in.(vars{j})=GetEndmember(fname,vars{j});
        end
    end
    
    for i=1:length(in.matdates)
        x=struct();
        [ declin, ~, solar_lon ]=Ephemeris(in.matdates(i)+10.5/24+8/24);
        mu0=sunang(lat,lon,declin,solar_lon);
%         mu=sunslope(mu0,phi0,Slope,Aspect);
        for j=1:length(vars)
            if j==4
            t=~isnan(x.grain_size) & x.grain_size>0;
            x.albedo=NaN(size(x.grain_size));    
                if dvisflag
                    x.albedo(t)=AlbedoLookup(double(x.grain_size(t)),...
                        double(mu0(t)),[],3,'dust',0);
                    x.albedo(t)=x.albedo(t)-0.63.*x.deltavis(t);
                else
                    x.albedo(t)=AlbedoLookup(double(x.grain_size(t)),...
                        double(mu0(t)),...
                        [],3,'dust',double(x.dust(t)).*1e-6);
                end
            else
            x.(vars{j})=rasterReprojection(in.(vars{j})(:,:,i),in.hdr.RefMatrix,...
                in.hdr.ProjectionStructure,target.ProjectionStructure,...
                'rasterref',target.RasterReference);
            end
            x.(vars{j})(isnan(x.(vars{j})))=0;
            x.(vars{j})(~mask)=NaN;
            axes(ha(j));
            imagesc(x.(vars{j}),'AlphaData',double(mask));
            text(1314,5,ltr{j},'FontSize',14,'VerticalAlignment','top',...
            'HorizontalAlignment','right');
            if j==1

            c1.Label.String=['fsca ' datestr(in.matdates(i))];
            end

        end
        frame=getframe(f1);
        writeVideo(f,frame)
    end
end
close(f);
close(f1);
