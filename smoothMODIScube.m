function [smoothedCube,refl,solarZ,cloudmask,pxweights]=...
    smoothMODIScube(tile,matdates,hdfdir,topofile,watermask)
%create time/space smoothed and gap filled (cloud-free) MOD09GA surface
%reflectance

%input:
%tile - tilename, e.g. 'h08v05'
%matdates - matdates for cube
%hdfdir - where the MOD09GA HDF files live for a certain tile, e.g. h08v04
%topofile- h5 file name from consolidateTopography, part of TopoHorizons
%watermask- logical mask w/ ones for water
%output:
%smoothedCube: smoothed and gap filled (cloud free) cube of MOD09GA values
%refl: terrain-corrected MOD09GA w/ NaNs for clouds and other missing data
%solarZ: solar zenith angles for cube
%cloudmask: cloudmask cube, logical
%snowmask: snowmask cube, logical
%pxweights: weight cube for each pixel (all bands together), 0-1

%get full directory listing for tile
d=dir(fullfile(hdfdir,['*.' tile '.*.hdf']));
d=struct2cell(d);
d=d(1,:);

[ ~,~,RasterReference] = sinusoidProjMODtile(tile);

nbands=7;
sz=[RasterReference.RasterReference_500m.RasterSize nbands length(matdates)];
refl=NaN(sz);
solarZ=NaN([sz(1) sz(2) sz(4)]);
cloudmask=false([sz(1) sz(2) sz(4)]);
pxweights=zeros([sz(1) sz(2) sz(4)]);
bandweights=zeros(sz);

% mask as cloud if swir band (6) is > than
swir_cloud_thresh=0.2;

parfor i=1:length(matdates)
    isodate=datenum2iso(matdates(i),7);
    m=regexp(d,['^MOD09GA.A' num2str(isodate) '\.*'],'once');
    m=~cellfun(@isempty,m);
    if any(m)
        f=fullfile(hdfdir,d{m});
        %get all band reflectance
        R=GetMOD09GA(f,'allbands'); 
        [~,aWeights,bWeights] = weightMOD09(f,topofile);
        pxweights(:,:,i)=aWeights;
        bandweights(:,:,:,i)=bWeights{1};
        
        x=single(GetMOD09GA(f,'SolarZenith'));
        if any(isnan(x(:)))
            x = inpaint_nans(double(x),4);
        end
        solarZ(:,:,i)=imresize(x,[sz(1) sz(2)]);

        x=single(GetMOD09GA(f,'SolarAzimuth'));
        if any(isnan(x(:)))
            x = inpaint_nans(double(x),4);
        end
        solarAzimuth=imresize(x,[sz(1) sz(2)]);

        %correct reflectance
        Rc=normalizeReflectance(R,topofile,solarZ(:,:,i),solarAzimuth);
        %fix negative values
        Rc(Rc<0.001)=0.001;
        %fix values > 1
        Rc(Rc>=1)=0.999;
        %create cloud & snow masks
        
        S=GetMOD09GA(f,'state');
        mod09gacm=imresize(S.cloud==1,[sz(1) sz(2)]);
        swir=Rc(:,:,6);       
        cm = mod09gacm & swir > swir_cloud_thresh;
        cloudmask(:,:,i)=cm;
        refl(:,:,:,i)=Rc;
        
        fprintf('loaded, corrected, and created masks for tile:%s date:%i \n',...
            tile,isodate);
    end
end

smoothedCube=NaN(size(refl));

for i=1:size(refl,3) % for each band
    tic;
    bandcube=squeeze(refl(:,:,i,:));
    bandcube(cloudmask)=NaN; %set all the clouds to NaN
    weights=squeeze(bandweights(:,:,i,:));
    weights(isnan(bandcube))=0;
    smoothedCube(:,:,i,:)=smoothDataCube(bandcube,weights,...
        'mask',~watermask,'method','smoothingspline');
    t2=toc;
    fprintf('filled and smoothed band:%i in %g min\n',i,t2/60);
end