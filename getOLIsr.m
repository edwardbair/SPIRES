function R=getOLIsr(ldir,target)
%retrieve OLI surface refl
%input: ldir - directory of SR tifs, string
%target - [] empty or target_hdr w/ fields RefMatrix and
%ProjectionStructure and rasterref
%output R - struct with fields bands, RefMatrix, ProjectionStructure,and
%RasterReference

%collection 1 on demand surface reflectance
d=dir(fullfile(ldir,'*band*.tif'));
% collection 1 ard surface reflectance
if isempty(d)
    d=dir(fullfile(ldir,'*SRB*.tif'));
end
if isempty(d)
    error('cannot find surface refl files');
end
for i=1:length(d)
    fname=fullfile(d(i).folder,d(i).name);
    X=single(geotiffread(fname));
    X(X==-9999)=NaN;
    X=X*1e-4;
    if i==1
        info=geotiffinfo(fname);
        RefMatrix=info.RefMatrix;
        ProjectionStructure=geotiff2mstruct(info);
        RasterReference=refmatToMapRasterReference(RefMatrix,size(X));
        if ~isempty(target)
            R.bands=zeros([target.RasterReference.RasterSize length(d)]);
        else
            R.bands=zeros([size(X(:,:,1)) length(d)]);
        end
    end
    
    if ~isempty(target)
        if any(RefMatrix(:)~=target.RefMatrix(:))
            % reproject if RefMatrices don't match
            [X,R.RefMatrix,R.RasterReference]=rasterReprojection(X,RefMatrix,...
                ProjectionStructure,target.ProjectionStructure,'rasterref',...
                target.RasterReference);
            R.ProjectionStructure=target.ProjectionStructure;
        else
            R.RefMatrix=RefMatrix;
            R.RasterReference=RasterReference;
            R.ProjectionStructure=ProjectionStructure;
        end
    end
    R.bands(:,:,i)=X;
end


end