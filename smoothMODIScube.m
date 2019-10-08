function [smoothedCube,refl,solarZ,cloudmask,pxweights]=...
    smoothMODIScube(tile,matdates,hdfdir,topofile,watermask)
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
% snowmask=false([sz(1) sz(2) sz(4)]);
cloudmask=false([sz(1) sz(2) sz(4)]);
pxweights=zeros([sz(1) sz(2) sz(4)]);
bandweights=zeros(sz);

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
%         [~,~,~,~, cloudScore, snowScore] =...
%         filterCloudSnow(snowR, Rc, solarZ(:,:,i));
        
        S=GetMOD09GA(f,'state');
        mod09gacm=imresize(S.cloud==1,[sz(1) sz(2)]);
        %cm =  cloudScore > 0 & mod09gacm;
%         red=Rc(:,:,1);
%         green=Rc(:,:,4);
%         blue=Rc(:,:,3);
        swir=Rc(:,:,6);
%         HOT=blue-0.5.*red-0.08;
%         NDSI=(green-swir)./(green+swir);
        
%         x=NDSI;
%         fillpx=isnan(x);
%         x(fillpx)=0;
%         o=runGabor(0:45:135,3,1,0.5,x);
%         o(fillpx)=0;
%         o2=rssq(o,3);
%         T=imadjust(o2,stretchlim(o2),[]);        
        cm = mod09gacm & swir > 0.2;
        cloudmask(:,:,i)=cm;
%         sm=snowScore >= 2 & ~cm;
%         snowmask(:,:,i)=sm;
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
    smoothedCube(:,:,i,:)=smoothDataCube(bandcube,weights,'mask',~watermask);
    t2=toc;
    fprintf('filled and smoothed band:%i in %g min\n',i,t2/60);
end