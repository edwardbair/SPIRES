function make_modis_video(infiles,target,subset)
%create reprojected MODIS video
%infiles - cell, list of h5 files to read
%output struct from smooth_and_run_modis
%w fields:
%fsca (canopy adj),
%grain radius (um)
%dust (ppmw)
%matdates (datenums)
%and hdr struct w it's own fields RefMatrix and ProjectionStructure
% target - hdr struct w. Projection Structure RasterReference to
%reproject to
%subset - matrix [row1 row2;col1 col2], e.g. [330 2400;1 2400]


fname='modis_fsca_modscag.avi';
f=VideoWriter(fname);

f.FrameRate=10;
f.Quality=90;
open(f);
vars={'snow_fraction','grain_size','dust'};

f1=figure('Position',[100 10 1500 750],'Color',[0.6 0.6 0.6]);
ha=tight_subplot(1, 3, 0.01, 0.01, 0);
set(ha,'NextPlot','replaceChildren');


for ii=1:length(infiles)
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
            axes(ha(j));
            ax=gca;
            ax.XAxis.Color = 'none';
            ax.YAxis.Color = 'none';
            set(ax,'XTick',[],'YTick',[],'YDir','reverse',...
                'Color',[0.6 0.6 0.6]);
            if j==1
                x(x==0)=NaN;
            end
            imagesc(x(subset(1,1):subset(1,2),subset(2,1):subset(2,2)));
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
end
close(f);
close(f1);
