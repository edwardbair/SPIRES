function make_modis_video(infiles,target,pshape)
%create reprojected MODIS video
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

fname='spires_video.avi';
f=VideoWriter(fname);

f.FrameRate=10;
f.Quality=90;
open(f);
vars={'snow_fraction','grain_size','dust'};

f1=figure('Position',[100 10 1500 750],'Color',[0.6 0.6 0.6]);
ha=tight_subplot(1, 3, 0.01, 0.01, 0);
set(ha,'NextPlot','replaceChildren');

[lon,lat]=pshape.boundary;

t=~isnan(lat) & ~isnan(lon);
[x,y]=mfwdtran(target.ProjectionStructure,lat,...
    lon);
[row,col]=map2pix(target.RasterReference,x,y);

mask=poly2mask(col(t),row(t),target.RasterReference.RasterSize(1),...
    target.RasterReference.RasterSize(2));

[lon,lat]=pshape.boundingbox;
[x,y]=mfwdtran(target.ProjectionStructure,lat,...
    lon);
[bbox_y,bbox_x]=map2pix(target.RasterReference,x,y);

% c=zeros(3,1);

cm=colormap(parula);

cm(1,:)=[0.4 0.4 0.4];

for j=1:3
    axes(ha(j));
    ax=gca;
    imagesc;
    colormap(cm);
%     plot(ha(j),col,row,'-k','LineWidth',0.25)
    xlim([bbox_x(1) bbox_x(2)+70]);
    ylim([bbox_y(2) bbox_y(1)+130]);
    
    ax.XAxis.Color = 'none';
    ax.YAxis.Color = 'none';
    set(ax,'XTick',[],'YTick',[],'YDir','reverse',...
        'Color',[0.6 0.6 0.6]);

     if j==1
        c1=colorbar('Location','south','Color','w');
        c1.Label.String='fsca';
        c1.Label.Color=[1 1 1];
        caxis([0 1]);
     elseif j==2
        c2=colorbar('Location','south','Color','w');
        c2.Label.String='grain radius, \mum';
        c2.Label.Color=[1 1 1];
        caxis([40 1200])
    elseif j==3
        c3=colorbar('Location','south','Color','w');
        c3.Label.String='dust conc, ppmw';
        c3.Label.Color=[1 1 1];
        caxis([0 500]);
    end
end

for ii=1:size(infiles,1)
    fname=infiles{ii};
    in=struct();
    for j=1:length(vars)
        if j==1
            [in.(vars{j}),in.matdates,in.hdr]=GetEndmember(fname,vars{j});
        else
            in.(vars{j})=GetEndmember(fname,vars{j});
        end
    end
    
    for i=1:length(in.matdates)
        for j=1:length(vars)
            x=rasterReprojection(in.(vars{j})(:,:,i),in.hdr.RefMatrix,...
                in.hdr.ProjectionStructure,target.ProjectionStructure,...
                'rasterref',target.RasterReference);
   
%             if j==1
% %                 x(x==0)=NaN;
%             end
%             x(~mask)=NaN;
            x(isnan(x))=0;
            x(~mask)=NaN;
            imagesc(ha(j),x,'AlphaData',double(mask));
            if j==1
%                 cc=colorbar('Location','south','Color','w');
        %        c.Label.String=['fsca ' datestr(in.matdates(i))];
            
            c1.Label.String=['fsca ' datestr(in.matdates(i))];
            end
%             freezeColors(ha(j),'nancolor',[0.6 0.6 0.6]);
        end
        frame=getframe(f1);
        writeVideo(f,frame)
    end
end
close(f);
close(f1);
