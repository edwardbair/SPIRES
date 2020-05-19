function animateMOD09GA(tiles,matdates,hdfbasedir,pshape,...
    target,bands)
%animate MOD09GA surface reflectance

%input:
%tiles - tilenames,  cell vector, e.g. {'h08v05','h08v04','h09v04'};
%matdates - matdates to process
%hdfbasedir - where the MOD09GA HDF files live must have sub directories
%pshape - boundary polyshape
%target - target hdr
%bands - vector to display, e.g. [1 4 3] for RGB

fname='modis_video.avi';
f=VideoWriter(fname);

f.FrameRate=10;
f.Quality=90;
open(f);

f1=figure('Position',[100 10 400 1000],'Color',[0.6 0.6 0.6]);
set(gca,'NextPlot','replaceChildren');

[lon,lat]=pshape.boundary;

t=~isnan(lat) & ~isnan(lon);
[x,y]=mfwdtran(target.ProjectionStructure,lat,lon);
[row,col]=map2pix(target.RasterReference,x,y);

mask=poly2mask(col(t),row(t),target.RasterReference.RasterSize(1),...
    target.RasterReference.RasterSize(2));
mask=repmat(mask,[1 1 length(bands)]);

[lon,lat]=pshape.boundingbox;
[x,y]=mfwdtran(target.ProjectionStructure,lat,lon);
[bbox_y,bbox_x]=map2pix(target.RasterReference,x,y);


ax=gca;
image;
xlim([bbox_x(1) bbox_x(2)+70]);
ylim([bbox_y(2) bbox_y(1)+130]);

ax.XAxis.Color = 'none';
ax.YAxis.Color = 'none';
set(ax,'XTick',[],'YTick',[],'YDir','reverse',...
    'Color',[0.6 0.6 0.6]);

if ~iscell(tiles)
    tiles={tiles};
end

% get all raster info
R=zeros([length(tiles),3,2]);
lr=zeros(length(tiles),2);
tsiz=zeros(length(tiles),2);

for i=1:length(tiles)
    [r,mstruct,rr] = sinusoidProjMODtile(tiles{i});
    R(i,:,:) = r.RefMatrix_500m;
    lr(i,:) = [rr.RasterReference_500m.RasterSize 1]*r.RefMatrix_500m;
    tsiz(i,:) = rr.RasterReference_500m.RasterSize;
end

%create refmat for mosaic
BigR=zeros(3,2);
BigR(3,1)=min(R(:,3,1));
BigR(2,1)=R(1,2,1); %assume same spacing for all pixels
BigR(1,2)=R(1,1,2);
BigR(3,2)=max(R(:,3,2));

%compute size of mosaic
%xy coords for lower rt corner
xy=[max(lr(:,1)) min(lr(:,2))];
sBig=map2pix(BigR,xy);
sBig=round(sBig);

sz=[sBig length(bands)];

%allocate cube
sz0=[target.RasterReference.RasterSize ...
    length(bands) length(matdates)];
refl=zeros(sz0);

parfor i=1:length(matdates)
    isodate=datenum2iso(matdates(i),7);
    %allocate daily cubes
    refl_=NaN(sz);
    %load up each tile
    for k=1:length(tiles)
        tile=tiles{k};
        
        %get full directory listing for tile
        d=dir(fullfile(hdfbasedir,tile,['*.' tile '.*.hdf']));
        d=struct2cell(d);
        d=d(1,:);
        assert(~isempty(d),'%s empty\n',hdfbasedir);
        m=regexp(d,['^MOD09GA.A' num2str(isodate) '\.*'],'once');
        m=~cellfun(@isempty,m);
        
        if any(m)
            [x,y]=pixcenters(squeeze(R(k,:,:)),tsiz(k,:));
            [r,c]=map2pix(BigR,x,y);
            r=round(r);
            c=round(c);
            f=fullfile(hdfbasedir,tile,d{m});
            %get all band reflectance
            for j=1:length(bands)
                sr=GetMOD09GA(f,['band' num2str(bands(j))]);
                refl_(r,c,j)=sr;
            end
            fprintf('loaded tile:%s date:%i \n',...
                tile,isodate);
        else
            fprintf('MOD09GA for %s %i not found,skipped\n',tile,isodate);
        end
    end
    refl_=rasterReprojection(refl_,BigR,mstruct,...
        target.ProjectionStructure,'rasterref',...
        target.RasterReference);
    refl_(~mask)=NaN;
    refl(:,:,:,i)=refl_;
end
for i=1:length(matdates)
    image(squeeze(refl(:,:,:,i)),'AlphaData',double(squeeze(mask(:,:,1))));
    %axis image;
    title(datestr(matdates(i)));
    frame=getframe(f1);
    writeVideo(f,frame)
end
close(f);
close(f1);
end
