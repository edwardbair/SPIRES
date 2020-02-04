function make_modis_video(in,target)
%create reprojected MODIS video
%in - output struct from smooth_and_run_modis
%w fields:
%fsca (canopy adj), 
%grain radius (um)
%dust (ppmw)
%matdates (datenums)
%and hdr struct w it's own fields RefMatrix and ProjectionStructure
% target - hdr struct w. Projection Structure RasterReference to
%reproject to
fname=[datestr(in.matdates(end),'yyyy') 'modis_fsca_modscag.avi'];
f=VideoWriter(fname);

f.FrameRate=10;
f.Quality=90;
open(f);
vars={'fsca','grainradius','dust'};

f1=figure('Position',[100 10 1500 750],'Color',[0.6 0.6 0.6]);
ha=tight_subplot(1, 3, 0.01, 0.01, 0);
set(ha,'NextPlot','replaceChildren');
for i=1:length(in.matdates)
    for j=1:length(vars)
        x=rasterReprojection(in.(vars{j})(:,:,i),in.hdr.RefMatrix,...
        in.hdr.ProjectionStructure,target.ProjectionStructure,...
        'rasterref',target.RasterReference);
        axes(ha(j));
        ax=gca;
        ax.XAxis.Color = 'none';
        ax.YAxis.Color = 'none';
        set(ax,'XTick',[],'YTick',[],'YDir','reverse',...
            'Color',[0.6 0.6 0.6]);
        if j==1
            x(x==0)=NaN;
        end
        imagesc(x(330:end,:));
        c=colorbar('Location','south','Color','w');
        freezeColors('nancolor',[0.6 0.6 0.6]);
        if j==1
            c.Label.String=['fsca ' datestr(in.matdates(i))];
            caxis([0 1]);
        elseif j==2
            c.Label.String='grain radius, \mum';
            caxis([0 1200])   
        elseif j==3
            c.Label.String='dust conc, ppmw';
            caxis([1 100]);
        end
        c.Label.Color=[1 1 1];
    end
    frame=getframe(f1);
    writeVideo(f,frame)
end
close(f);
close(f1);
