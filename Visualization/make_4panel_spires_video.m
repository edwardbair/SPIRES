function make_4panel_spires_video(infiles,vidname,rrb)
%create reprojected spires MODIS video
%input: infile - input h5 files, cell Nx1
%vidname : output vidname
%rrb : target rasterref w/ CRS
f=VideoWriter(vidname);

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

f1=figure('Position',[100   1   1000   1000],'Color','w');
set(f1,'toolbar','none');
ha=tight_subplot(2, 2, [0.025 0.001], [0.025 0.01], [0.025 0]);
set(ha,'NextPlot','replaceChildren');

cm=colormap(parula);

cm(1,:)=[0.5 0.5 0.5];

for j=1:length(vars)
    axes(ha(j));
    ax=gca;
    imagesc;
    axis image;
    
    colormap(cm);
    
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
%     set(gca,'YDir','reverse','Box','on','XTicklabel',[],'YTickLabel',[]);
%     set(gca,'YTick',rlat,'YTickLabel',num2str(lat_l'));
%     set(gca,'XTick',clon,'XTickLabel',num2str(lon_l'))
    
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
        caxis([0 1]);
     elseif j==2
        c2=colorbar('Location','south');
        c2.Label.String='grain radius, \mum';
        caxis([40 900])
    elseif j==3
        c3=colorbar('Location','south');
        if dvisflag
            c3.Label.String='deltavis';
            caxis([0 0.4]);
        else
            c3.Label.String='dust conc, ppm';
            caxis([0 300]);
        end
     elseif j==4
        c4=colorbar('Location','south');
        c4.Label.String='albedo';
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
        [ declin, ~, solar_lon ]=EarthEphemeris(in.matdates(i)+10.5/24+8/24);
        mu0=sunang(lat,lon,declin,solar_lon);

        for j=1:length(vars)
            if j==4
            t=~isnan(x.grain_size) & x.grain_size>0;
            x.albedo=NaN(size(x.grain_size));    
                if dvisflag
                    x.albedo(t)=AlbedoLookup(double(x.grain_size(t)),...
                        double(mu0(t)),3,'dust',0);
                    x.albedo(t)=x.albedo(t)-0.63.*x.deltavis(t);
                else
                    x.albedo(t)=AlbedoLookup(double(x.grain_size(t)),...
                        double(mu0(t)),...
                        3,'dust',double(x.dust(t)).*1e-6);
                end
            else
            x.(vars{j})=rasterReprojection(in.(vars{j})(:,:,i),...
                in.hdr.RasterReference,'InProj',in.hdr.ProjectionStructure,...
                'rasterref',rrb);
            x.(vars{j})(x.(vars{j})<0)=0;
            
            end
%            
%             x.(vars{j})(~mask)=NaN;
            axes(ha(j));
            if j==1
                mask=~isnan(x.(vars{1}));
            else
                mask=~isnan(x.(vars{1}));
                x.(vars{j})(isnan(x.(vars{j})))=0;
            end
            imagesc(x.(vars{j}),'AlphaData',mask);
            text(1,1,ltr{j},'FontSize',14,'VerticalAlignment','top',...
            'HorizontalAlignment','right','Units','normalized');
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
