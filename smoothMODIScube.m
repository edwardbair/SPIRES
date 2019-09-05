function [smoothedCube,refl,solarZ]=smoothMODIScube(tile,matdates,...
    hdfdir,topofile,snowR)
%create time/space smoothed and gap filled (cloud-free) MOD09GA surface
%reflectance
%input:
%tile - tilename, e.g. 'h08v05'
%matdates - matdates for cube
%hdfdir - where the MOD09GA HDF files live for a certain tile, e.g. h08v04
%topofile- h5 file name from consolidateTopography, part of TopoHorizons
%snowR - snow reflectane structure created using prepCloudSnowFilter
%output:
%smoothedCube: smoothed and gap filled (cloud free) cube of MOD09GA values
%refl: terrain-corrected MOD09GA w/ NaNs for clouds and other missing data
%solarZ: solar zenith angles for cube

%get full directory listing for tile
d=dir(fullfile(hdfdir,['*.' tile '.*.hdf']));
d=struct2cell(d);
d=d(1,:);

[ RefMatrix,ProjectionStructure,RasterReference] = sinusoidProjMODtile(tile);

nbands=7;
sz=[RasterReference.RasterReference_500m.RasterSize nbands length(matdates)];
refl=NaN(sz);
solarZ=NaN([sz(1) sz(2) sz(4)]);
%snowmask=false([sz(1) sz(2) sz(4)]);
cloudmask=false([sz(1) sz(2) sz(4)]);
bandweights=zeros(sz);

parfor i=1:length(matdates)
    isodate=datenum2iso(matdates(i),7);
    m=regexp(d,['^MOD09GA.A' num2str(isodate) '\.*'],'once');
    m=~cellfun(@isempty,m);
    if any(m)   
        f=fullfile(hdfdir,d{m});
        %get 7 band reflectance
        R=GetMOD09GA(f,'allbands'); 
        [~,~,bWeights] = weightMOD09(f,topofile);
        bandweights(:,:,:,i)=bWeights{1};
        
        %get illumination data and
%         x=GetMOD09GA(f,'SensorZenith');
%         %reproject to 500 m
%         sensorZ(:,:,i)=rasterReprojection(x,RefMatrix.RefMatrix_1km,...
%             ProjectionStructure,ProjectionStructure,'rasterref',...
%             RasterReference.RasterReference_500m);
        x=single(GetMOD09GA(f,'SolarZenith'));
        solarZ(:,:,i)=rasterReprojection(x,RefMatrix.RefMatrix_1km,...
            ProjectionStructure,ProjectionStructure,'rasterref',...
            RasterReference.RasterReference_500m);
        x=single(GetMOD09GA(f,'SolarAzimuth'));
        solarAzimuth=rasterReprojection(x,RefMatrix.RefMatrix_1km,...
            ProjectionStructure,ProjectionStructure,'rasterref',...
            RasterReference.RasterReference_500m);
        %correct reflectance
        Rc=normalizeReflectance(R,topofile,solarZ(:,:,i),solarAzimuth);
        %put NaNs in for clouds
        [likelySnow, maybeSnow, likelyCloud, maybeCloud ] =...
        filterCloudSnow(snowR, Rc, solarZ(:,:,i));
        %set likelyCloud and MaybeCloud to -1
%         sm=likelySnow | maybeSnow;
%         snowmask(:,:,i)=sm;
        cm=likelyCloud | maybeCloud;
        cloudmask(:,:,i)=cm;
        refl(:,:,:,i)=Rc;
    end
end
   
%fix negative values
refl(refl<0.001)=0.001;
toomanycloudydays=7;
cloudmask=cloudPersistenceFilter(cloudmask,toomanycloudydays);
% toofewsnowdays=3;
% snowmask=snowPersistenceFilter(snowmask,toofewsnowdays,1);

smoothedCube=NaN(size(refl));

for i=1:size(refl,3)
    bandcube=squeeze(refl(:,:,i,:));
    bandcube(cloudmask)=NaN; %set all the clouds to NaN
    %bandcube=infillDataCube(bandcube);
    weights=squeeze(bandweights(:,:,i,:));
    weights(isnan(bandcube))=0;
    smoothedCube(:,:,i,:)=smoothDataCube(bandcube,weights);
end