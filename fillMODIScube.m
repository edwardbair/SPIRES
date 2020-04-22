function [filledCube,refl,solarZ,cloudmask,pxweights]=...
    fillMODIScube(tile,matdates,hdfdir,topofile,watermask)
%create gap filled (e.g. cloud-free) MOD09GA surface
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
% bandweights=zeros(sz);

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
        [~,aWeights,~] = weightMOD09(f,topofile);
        pxweights(:,:,i)=aWeights;
%         bandweights(:,:,:,i)=bWeights{1};
        
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
        Slope=GetTopography(topofile,'slope');
        Aspect=GetTopography(topofile,'aspect');
        Rc=normalizeReflectance(R,Slope,Aspect,solarZ(:,:,i),solarAzimuth);
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

filledCube=NaN(size(refl));
wm=reshape(watermask,[sz(1)*sz(2) 1]);

for i=1:size(refl,3) % for each band
    tic;
    bandcube=squeeze(refl(:,:,i,:));
    bandcube(cloudmask)=NaN; %set all the clouds to NaN
    %fill datacube
    vec=reshape(bandcube,[sz(1)*sz(2) sz(4)])'; %col major (days x pixels)
    parfor j=1:size(vec,2)
        if ~wm(j)
            v=vec(:,j)
            t=isnan(v);
            if nnz(~t) > 3
                %interpolate
                v(t)=interp1(matdates(~t),v(~t),matdates(t));
                %extrapolate
                tt=isnan(v);
                if any(tt)
                    v(tt)=interp1(matdates(~tt),v(~tt),matdates(tt),...
                        'nearest','extrap');
                end
                vec(:,j)=v;
            end
        end
    end
    XX=reshape(vec',[sz(1) sz(2) sz(4)]);
%     t=isnan(XX);
%     XX(t)=0; %set all remaining NaNs to zero
    parfor j=1:size(XX,3)
        if any(isnan(XX(:,:,j)) & ~watermask,'all') %fill any remaining NaNs 
            %that could  not be temporally interpolated spatially
             XX(:,:,j)=inpaint_nans(double(XX(:,:,j)),4);
        end
    end
    filledCube(:,:,i,:)=XX;
    t2=toc;
    fprintf('filled band:%i in %g min\n',i,t2/60);
end