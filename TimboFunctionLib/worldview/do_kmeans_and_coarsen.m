%paths
Kings_dir='C:\raid\data\nbair\worldview\20150531_r3';
outdir=fullfile(Kings_dir,'3m');
d=dir(fullfile(Kings_dir,'*R*C*.tif'));
asodir='C:\raid\data\nbair\datasets\ASO\data\kings_basinCA\snowon\USCAKC20150531';
target_info=geotiffinfo(fullfile(asodir,'KC20150531_SUPERsnow_depth_MANUAL.tif'));
target_mstruct=geotiff2mstruct(target_info);
addpath(genpath('C:\raid\data\jeff\OneDrive\MATLAB\FunctionLib'));
addpath('C:\raid\data\nbair\MATLAB functions');
%% reproject to intermediate pixel size geog.
parfor i=1:length(d);
    fname=fullfile(Kings_dir,d(i).name);
    tic
    fprintf('starting %s\n',fname)
    x=geotiffread(fname);
    x=single(x);
    ind=x(:,:,:)==0;
    x(ind)=NaN;
    info=geotiffinfo(fname);    
    RefMatrix=info.RefMatrix;
    rr=refmatToGeoRasterReference(RefMatrix,[info.Height info.Width]);
    %use a pixel size that is smaller than target and an even number
    [x2,r2,rr2]=reprojectRaster(x,RefMatrix,[],[],'pixelsize',...
        3e-6*6,'adjust',true,'method','bicubic');
    savename=fullfile(outdir,strcat('3m_',d(i).name));
    geotiffwrite(savename,uint8(x2),rr2,'GeoKeyDirectoryTag',...
        info.GeoTIFFTags.GeoKeyDirectoryTag);
    fprintf('finished %s\n',fname)
    toc
end
%% mosaic
d=dir(fullfile(outdir,'3m_*.tif'));
filenames=cell(length(d),1);
for j=1:length(d);
    filenames{j}=fullfile(outdir,d(j).name);
end
[Big,BigRmap]=mosaicTiles(filenames,'tif','uint8');
save([outdir,'\','Big.mat'],'Big','BigRmap');
%% fill in seams and reproject to UTM
[X,Y]=meshgrid(1:size(Big,2),1:size(Big,1));
for i=1:size(Big,3);
    bigchannel=Big(:,:,i);
    ind=bigchannel==0;
    se=strel('square',10);
    ind2=imdilate(ind,se);
    Vq=scatteredInterpolant(X(ind2 & ~ind),Y(ind2 & ~ind),...
        double(bigchannel(ind2 & ~ind)),'linear','none');
    bigchannel(ind)=Vq(X(ind),Y(ind));
    Big(:,:,i)=bigchannel;
end
save([outdir,'BigInterp.mat'],'Big','BigRmap');
%reproject to UTM
rr=refmatToMapRasterReference(target_info.RefMatrix,...
    [target_info.Height target_info.Width]);
[Big_re,Big_re_Refmatrix]=reprojectRaster(Big,BigRmap,...
    [],target_mstruct,'rasterref',rr);
save([outdir,'\','Big_re.mat'],'Big_re','Big_re_Refmatrix');
%% cut out tile shape & Kings basin shape
%tile shape
c=shaperead(fullfile(Kings_dir,'20150531c_SEAMLINES_SHAPE',...
    '20150531c_SEAMLINES_SHAPE.shp'));
b=shaperead(fullfile(Kings_dir,'20150531b_SEAMLINES_SHAPE',...
    '20150531b_SEAMLINES_SHAPE.shp'));
a=shaperead(fullfile(Kings_dir,'20150531a_SEAMLINES_SHAPE',...
    '20150531a_SEAMLINES_SHAPE.shp'));

[x,y]=polybool('union',c.X,c.Y,b.X,b.Y);
[x,y]=polybool('union',x,y,a.X,a.Y);

[x,y]=projfwd(target_mstruct,y,x);
[row,col]=map2pix(Big_re_Refmatrix,x,y);
row=round(row);
col=round(col);

seam_mask=poly2mask(col,row,sBig(1),sBig(2));
seam_mask=repmat(seam_mask,[1 1 3]);
Big_re(~seam_mask)=0;

%basin shape
KB=shaperead('C:\raid\data\nbair\Basins\Kings20150403.shp');

[x,y]=projfwd(target_mstruct,KB.Y,KB.X);
[row,col]=map2pix(Big_re_Refmatrix,x,y);
row=round(row);
col=round(col);
t=isnan(row) | isnan(col);
row(t)=[];
col(t)=[];

basin_mask=poly2mask(col,row,sBig(1),sBig(2));
basin_mask=repmat(basin_mask,[1 1 3]);

%set mosaic to zero outside basin boundary 
Big_re(~basin_mask)=0;
geotiffwrite([outdir,'\','mosaic_3m.tif'],Big_re,Big_re_Refmatrix,...
    'GeoKeyDirectoryTag',target_info.GeoTIFFTags.GeoKeyDirectoryTag);
%% kmeans
sBig=size(Big_re);
Big_re=single(Big_re);
ind=Big_re(:,:,:)==0;

null_vec=reshape(ind,[sBig(1)*sBig(2) sBig(3)]);
null_vec=~any(~null_vec,2);

Big_vec=reshape(Big_re,[sBig(1)*sBig(2) sBig(3)]);
Big_vec(null_vec,:)=[];

idx=kmeans(Big_vec,5,'MaxIter',500,'Display','iter',...
    'OnlinePhase','off');
idxfull=zeros([sBig(1)*sBig(2) 1],'uint8');
idxfull(~null_vec)=idx;
idxfull=reshape(uint8(idxfull),sBig(1),sBig(2));
save([outdir,'\','kmeans_3m.mat'],'idxfull');
geotiffwrite([outdir,'\','kmeans_3m.tif'],idxfull,Big_re_Refmatrix,...
    'GeoKeyDirectoryTag',target_info.GeoTIFFTags.GeoKeyDirectoryTag);
