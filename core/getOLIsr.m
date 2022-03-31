function [R,cosm]=getOLIsr(ldir,target)
%retrieve OLI surface refl
%deals with collection 1 (on demand or ard),2, or HLS (S30)
%input: ldir - directory of SR tifs, string
%target - [] empty or target_hdr w/ fields RefMatrix and
%ProjectionStructure and rasterref
%output:
%   R - struct with fields bands, RefMatrix, ProjectionStructure,and
%   RasterReference
%   cosm - cloud or snow mask from BQA data.(use later to remove bright,
%   warm, non snow pixels)

%directory listing produces ascending sort of bands, except for HLS
%collection 1 on demand surface reflectance
d=dir(fullfile(ldir,'*band*.tif'));
product='collection1';
% collection 1 ard surface reflectance
if isempty(d)
    d=dir(fullfile(ldir,'*SRB*.tif'));
    product='collection1';
end
%collection 2 surface reflectance
if isempty(d)
    d=dir(fullfile(ldir,'*_SR_B*.TIF'));
    product='collection2';
end
%HLS S30 surface reflectance
if isempty(d)
    d=dir(fullfile(ldir,'HLS*v*B*.tif'));
    idx=[1:8 13 9:12]; %puts bands in ascending order of bandpass and discards angle files
    d=d(idx);
    product='HLS';
end
if isempty(d)
    error('cannot find surface refl files');
end
for i=1:length(d)
    fname=fullfile(d(i).folder,d(i).name);
    X=single(readgeoraster(fname));

    switch product
        case {'collection1','HLS'}
            X(X==-9999)=NaN;
            X=X*1e-4;
        case 'collection2'
            if i==1
                bqafname=d(i).name;
                bqaname=strrep(bqafname,'SR_B1','QA_PIXEL');
                bqa=readgeoraster(fullfile(d(i).folder,bqaname));
                S=unpackLandsat8BQA(bqa,'collection2');
            end
            X(S.fill)=NaN;
            X=X*2.75e-5-0.2; %rescale
    end
    if i==1
        %         switch product
        %             case 'collection1'
        %                 info=geotiffinfo(fname);
        %                 ProjectionStructure=geotiff2mstruct(info);
        %                 RefMatrix=info.RefMatrix;
        %             case {'collection2','HLS'} %build mstruct from CRS info
        info=georasterinfo(fname);
        RasterReference=info.RasterReference;
        %                 ProjectionStructure=defaultm('tranmerc');
        %                 ProjectionStructure.falseeasting=...
        %                     info.CoordinateReferenceSystem.ProjectionParameters.FalseEasting;
        %                 ProjectionStructure.falsenorthing=...
        %                     info.CoordinateReferenceSystem.ProjectionParameters.FalseNorthing;
        %                 ProjectionStructure.geoid=...
        %                     [info.CoordinateReferenceSystem.GeographicCRS.Spheroid.SemimajorAxis ...
        %                     info.CoordinateReferenceSystem.GeographicCRS.Spheroid.Eccentricity];
        %                 ProjectionStructure.origin=[...
        %                     info.CoordinateReferenceSystem.ProjectionParameters.LatitudeOfNaturalOrigin ...
        %                     info.CoordinateReferenceSystem.ProjectionParameters.LongitudeOfNaturalOrigin ...
        %                     0];
        %                 ProjectionStructure.scalefactor=...
        %                     info.CoordinateReferenceSystem.ProjectionParameters.ScaleFactorAtNaturalOrigin;
        %                 RefMatrix=RasterRef2RefMat(info.RasterReference);
        %         end
        %         RasterReference=refmatToMapRasterReference(RefMatrix,size(X));
        if ~isempty(target)
            R.bands=zeros([target.RasterReference.RasterSize length(d)]);
        else
            R.bands=zeros([size(X(:,:,1)) length(d)]);
        end
    end

    if ~isempty(target)
        %         if any(RefMatrix(:)~=target.RefMatrix(:)) || ...
        if any(size(X)~=target.RasterReference.RasterSize)

            % reproject if RefMatrices or raster sizes don't match
            %             [X,R.RasterReference]=rasterReprojection(X,RasterReference,...
            %                 'InProj',ProjectionStructure,'OutProj',target.ProjectionStructure,'rasterref',...
            %                 target.RasterReference);
            %             R.RefMatrix=RasterRef2RefMat(R.RasterReference);
            %             R.ProjectionStructure=target.ProjectionStructure;
            X=rasterReprojection(X,RasterReference,...
                'rasterref',target.RasterReference,'cells',false);
        end
%     else
        %         R.RefMatrix=RefMatrix;
%         R.RasterReference=RasterReference;
        %         R.ProjectionStructure=ProjectionStructure;
    end
    R.bands(:,:,i)=X;
end

%output BQA mask for fsca adjustment, reproject if needed
switch product
    case {'collection1','HLS'}
        %not yet written
        cosm=[];
    case 'collection2'
        %S structure contains BQA info
        %cloud or snow mask
        %conservative cloud mask
        cm=S.cloudConfidence==2|S.cloudConfidence==3;
        cm=cm|S.cirrus|S.cloud|S.dilatedCloud;
        sm=S.snow;

        %buffer snow mask for low fsca around edges of snowfields (~1k buffer)
        SE1=strel('disk',35);
        %SE2=strel('disk',5);

        sm=imdilate(sm,SE1);
        cosm= cm | sm; %& ~S.water;mask water later
        %fill small holes in the cosm mask (~600m holes)
        cosm = ~bwareaopen(~cosm, 20);

        if ~isempty(target)
            %             if any(RefMatrix(:)~=target.RefMatrix(:)) || ...
            if any(size(cosm)~=target.RasterReference.RasterSize)
                % reproject if RefMatrices or raster sizes don't match
                cosm=rasterReprojection(cosm,RasterReference,...
                    'rasterref',target.RasterReference,'method','nearest','cells',false);
            end
        end
end
R.RasterReference=RasterReference;
end