function make_4panel_spires_video_p(infiles,vidname,rrb,varargin)
%parallel implementation of function to create reprojected spires MODIS video
%input: infile - input h5 files, cell Nx1
%vidname : output vidname
%rrb : target rasterref w/ CRS
%optional: mask matching rrb, logical

mask0=[];
if nargin==4
    mask0=varargin{1};
end

vars={'snow_fraction','grain_size','dust','albedo'};
info=h5info(infiles{1});

dvisflag=false;
for i=1:length(info.Groups.Groups.Datasets)
    if strcmp(info.Groups.Groups.Datasets(i).Name,'deltavis')
        dvisflag=true;
        vars={'snow_fraction','grain_size','deltavis','albedo'};
    end
end

[x,y]=worldGrid(rrb);
[lat,lon]=projinv(rrb.ProjectedCRS,x,y);

%     dy=-2;
%     dx=2;
%     y1=ceil(lat(1,1));
%     y2=floor(lat(end,1));
%     x1=ceil(lon(1,1));
%     x2=floor(lon(1,end));    
%     lat_l=y1:dy:y2;
%     lon_l=x1:dx:x2;

    y1=lat(1,1);
    y2=lat(end,1);
    x1=lon(1,1);
    x2=lon(1,end);    
    lat_l=linspace(y1,y2,4);
    lon_l=linspace(x1,x2,4);
    
    [x,y]=projfwd(rrb.ProjectedCRS,mean(lat_l)*ones(size(lon_l)),lon_l);
    [clon,~]=worldToIntrinsic(rrb,x,y);
    
    [x,y]=projfwd(rrb.ProjectedCRS,lat_l,mean(lon_l)*ones(size(lat_l)));
    [~,rlat]=worldToIntrinsic(rrb,x,y);

cm=colormap(parula);
close;
cm(1,:)=[0.5 0.5 0.5];
%ltr={'(a)','(b)','(c)','(d)'};

spmd
    figure('Position',[1   1   1900 1200],'Color','w','Visible','off');
    set(gcf,'toolbar','none');
    tiledlayout(2,2,'TileSpacing','none','padding','tight');
    for j=1:length(vars)
        nexttile(j)
        imagesc;
        axis image;
        colormap(cm);
        set(gca,'YDir','reverse','Box','on','XTick',clon,'XTicklabel',[],...
            'YTick',rlat,'YTickLabel',[],'box','off');
        set(gca,'NextPlot','replaceChildren');
        if j==1 || j==3
            set(gca,'YTick',rlat,'YTickLabel',num2str(lat_l'));
        end
        if j>2
            set(gca,'XTick',clon,'XTickLabel',num2str(lon_l'))
        end
        c=colorbar('Location','south');
%         c.Position(1)=c.Position(1)+0.25;
        c.Position(2)=c.Position(2)+0.05;
        c.Position(3)=c.Position(3)-0.16;
        if j==1
            %colorbar('Location','south');
            c.Label.String='fsca';
            caxis([0 1]);
        elseif j==2
            %colorbar('Location','south');
            c.Label.String='grain radius, \mum';
            caxis([40 900])
        elseif j==3
            %colorbar('Location','south');
            if dvisflag
                c.Label.String='deltavis';
                caxis([0 0.4]);
            else
                c.Label.String='dust conc, ppm';
                caxis([0 300]);
            end
        elseif j==4
            %colorbar('Location','south');
            c.Label.String='albedo';
            caxis([0.4 0.9]);
        end
    end
    c.FontSize=25;
end
frames={};

for ii=1:size(infiles,1)
    fname=infiles{ii};
    for j=1:3
        if j==1
            [fsca,matdates,hdr]=GetEndmember(fname,vars{j});
        elseif j==2
            grain_size=GetEndmember(fname,vars{j});
        elseif j==3
            dust=nan([1 1 length(matdates)]);
            deltavis=nan([1 1 length(matdates)]);
            if dvisflag
                deltavis=GetEndmember(fname,vars{j});
            else
                dust=GetEndmember(fname,vars{j});
            end
        end
    end

    parfor i=1:length(matdates)
        %need to initialize as empty temp var for parfor
        x=[];
        fscaf=[];
        grain_sizef=[];
        dustf=[];
        deltavisf=[];
        mask=[];

        fsca_i=fsca(:,:,i);

        grain_size_i=grain_size(:,:,i);
        deltavis_i=[];
        dust_i=[];
        if dvisflag
            deltavis_i=deltavis(:,:,i);
        else
            dust_i=dust(:,:,i);
        end

        matdates_i=matdates(i);

        for j=1:length(vars)
            if j==1
                fscaf=rasterReprojection(fsca_i,...
                    hdr.RasterReference,'InProj',hdr.ProjectionStructure,...
                    'rasterref',rrb);
                fscaf(fscaf<0)=0;
                if ~isempty(mask0)
                    mask=mask0 & ~isnan(fscaf);
                else
                    mask=~isnan(fscaf);
                end
                x=fscaf;
            elseif j==2
                grain_sizef=rasterReprojection(grain_size_i,...
                    hdr.RasterReference,'InProj',hdr.ProjectionStructure,...
                    'rasterref',rrb);
                grain_sizef(grain_sizef<0)=0;
                grain_sizef(isnan(grain_sizef))=0;
                x=grain_sizef;
            elseif j==3
                if dvisflag
                    deltavisf=rasterReprojection(deltavis_i,...
                        hdr.RasterReference,'InProj',hdr.ProjectionStructure,...
                        'rasterref',rrb);
                    deltavisf(deltavisf<0)=0;
                    deltavisf(isnan(deltavisf))=0;
                    x=deltavisf;
                else
                    dustf=rasterReprojection(dust_i,...
                        hdr.RasterReference,'InProj',hdr.ProjectionStructure,...
                        'rasterref',rrb);
                    dustf(dustf<0)=0;
                    dustf(isnan(dustf))=0;
                    x=dustf;
                end
            elseif j==4
                [ declin, ~, solar_lon ]=EarthEphemeris(matdates_i+10.5/24+8/24);
                mu0=sunang(lat,lon,declin,solar_lon);

                t=~isnan(grain_sizef) & grain_sizef>0;
                
                albedo=NaN(size(grain_sizef));
                if any(t,'all')
                    if dvisflag
                        albedo(t)=AlbedoLookup(double(grain_sizef(t)),...
                            double(mu0(t)),3,'dust',0);
                        albedo(t)=albedo(t)-0.63.*deltavisf(t);
                    else
                        albedo(t)=AlbedoLookup(double(grain_sizef(t)),...
                            double(mu0(t)),...
                            3,'dust',double(dustf(t)).*1e-6);
                    end
                end
                x=albedo;

            end
            nexttile(j)
            imagesc(x,'AlphaData',mask);

%             text(1,0.5,ltr{j},'FontSize',25,'VerticalAlignment','bottom',...
%                 'HorizontalAlignment','right','Units','normalized');
            if j==1
                text(0,0.25,datestr(matdates_i),'units','normalized',...
                    'FontSize',25);
            end
        set(gca,'FontSize',25);
        end
%         child=get(gcf,'children');
%         child=child.Children;
%         for y = 1:length(child)
%             chi=child(y);
%             set(chi, 'fontsize', 25);
%         end


        frame=getframe(gcf);
        %parfor ensures frame order will be correct
        frames= [frames, frame];
        fprintf('done w/ %s\n', datestr(matdates_i));
    end
end

f=VideoWriter(vidname);
f.FrameRate=10;
f.Quality=90;
open(f);

for idx=1:numel(frames)
    writeVideo(f,frames{idx})
end

close(f);