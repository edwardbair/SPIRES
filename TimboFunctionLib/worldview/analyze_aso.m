dir_3m='C:\raid\data\nbair\datasets\worldview\20150531_r3\3m';
aso_dir='C:\raid\scratch\nbair\ASO\data\kings_basinCA\snowon\USCAKC20150531';
km=load(fullfile(dir_3m,'kmeans_3m.mat'));
classes=km.idxfull;
aso=geotiffread(fullfile(aso_dir,'KC20150531_SUPERsnow_depth_MANUAL.tif'));
info=geotiffinfo(fullfile(aso_dir,'KC20150531_SUPERsnow_depth_MANUAL.tif'));
mstruct=geotiff2mstruct(info);
kingsShape=shaperead(...
    'C:\raid\data\nbair\datasets\Basins\Kings20150403.shp');
[x,y]=projfwd(mstruct,kingsShape.Y,kingsShape.X);
[row,col]=map2pix(info.RefMatrix,x,y);
t=isnan(row) | isnan(col);
row(t)=[];
col(t)=[];
mask=poly2mask(col,row,info.Height,info.Width);

nodata= aso < 0 | classes == 0 | ~mask;
aso=aso>0;
snow=classes==2;
snow(nodata)=false;
aso(nodata)=false;

true_positive = aso & snow & ~nodata;
false_positive = aso & ~snow & ~nodata;
true_negative= ~aso & ~snow & ~nodata;
false_negative = ~aso & snow & ~nodata;

%take best of 3x3 neighborhood
true_positive=ordfilt2(true_positive,9,ones(3,3));
false_positive=ordfilt2(false_positive,1,ones(3,3));
true_negative=ordfilt2(true_negative,9,ones(3,3));
false_negative=ordfilt2(false_negative,1,ones(3,3));
    
c_table=array2table(...
    [nnz(true_positive) nnz(false_positive); nnz(false_negative) nnz(true_negative)],...
    'RowNames',{'ASO_snow';'ASO_no_snow'},'VariableNames',...
    {'WV_snow','WV_no_snow'});
N=sum(c_table{:,1})+sum(c_table{:,2});
Recall=nnz(true_positive)/(nnz(true_positive)+nnz(false_positive)); % also frequency of hits
Precision=nnz(true_positive)/(nnz(true_positive)+nnz(false_negative)); % also prob. detection/sensitivity

Accuracy=(nnz(true_positive)+nnz(true_negative))/...
    (nnz(true_positive)+nnz(true_negative)+nnz(false_positive)+nnz(false_negative));
Fstat=2*(Recall*Precision)/(Recall+Precision);

%Worldview-3 VNIR, Accuracy: RMSE 3.9 m, CE-90 accuracy 8.9 m

geotiffwrite(fullfile(dir_3m,'fp.tif'),false_positive,info.RefMatrix,'GeoKeyDirectoryTag',...
    info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(fullfile(dir_3m,'fn.tif'),false_negative,info.RefMatrix,'GeoKeyDirectoryTag',...
    info.GeoTIFFTags.GeoKeyDirectoryTag);

