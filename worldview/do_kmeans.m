Kings_dir='C:\raid\data\nbair\worldview\20150531Kings';
d=dir(fullfile(Kings_dir,'*R*C*.tif'));

k=3; %classes
parfor i=1:length(d);
    fname=fullfile(Kings_dir,d(i).name);
    tic
    fprintf('starting %s\n',fname)
    x=geotiffread(fname);
    info=geotiffinfo(fname);
    [nrows,ncols,C]=size(x);
    RefMatrix=info.RefMatrix;
    rasterref=refmatToGeoRasterReference(RefMatrix,[nrows ncols]);
    x=reshape(x,nrows*ncols,C);
    idx=kmeans(double(x),4,'Display','iter');
    idx=reshape(idx,nrows,ncols);
    savename=fullfile(Kings_dir,strcat('kmeans',d(i).name));
    geotiffwrite(savename,idx,rasterref,'GeoKeyDirectoryTag',...
        info.GeoTIFFTags.GeoKeyDirectoryTag);
    fprintf('finished %s\n',fname)
    toc
end