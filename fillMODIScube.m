function [filledCube,refl,SolarZenith,cloudmask,pxweights]=...
    fillMODIScube(tiles,matdates,hdfbasedir,topodir,topofile,mask,swir_b)
%create gap filled (e.g. cloud-free) MOD09GA surface
%reflectance

%input:
%tiles - tilenames,  cell vector, e.g. {'h08v05','h08v04','h09v04'};
%matdates - matdates to process
%hdfbasedir - where the MOD09GA HDF files live must have sub directories 
% that correspond to entries in tile, e.g. 'h08v04'
%topodir - directory for h5 topo files from from consolidateTopography, 
%part of TopoHorizon that contain topofiles for each tile in tiles, 
%e.g. h09v05dem_463m_Topography.h5
%topofile - h5 target topo file name. This topography file contains the target
%geographic information that everything will be tiled and
%cropped/reprojected to
%mask - logical mask w/ ones for areas to exclude, same geog extent as
%topofile
%swir_b - swir band, scalar
%output:
%filledCube: gap filled (cloud free) cube of MOD09GA values
%refl: terrain-corrected MOD09GA w/ NaNs for clouds and other missing data
%solarZ: solar zenith angles for cube
%cloudmask: cloudmask cube, logical
%pxweights: weight cube for each pixel (all bands together), 0-1

nbands=7;
% mask as cloud if swir band (6) is > than
swir_cloud_thresh=0.2;

if ~iscell(tiles)
    tiles={tiles};
end

% get all raster info
R=zeros([length(tiles),3,2]);
lr=zeros(length(tiles),2);
tsiz=zeros(length(tiles),2);

for i=1:length(tiles)
    [r,mstruct,rr] = sinusoidProjMODtile(tiles{i});
    R(i,:,:)=r.RefMatrix_500m;
    lr(i,:)=[rr.RasterReference_500m.RasterSize 1]*r.RefMatrix_500m;
    tsiz(i,:)=rr.RasterReference_500m.RasterSize;
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

sz=[sBig nbands];
% sunnames={'SolarZenith','SolarAzimuth'};

[Slope,hdr]=GetTopography(topofile,'slope');
Aspect=GetTopography(topofile,'aspect');

% [cropr,cropc]=map2pix(BigR,hdr.RasterReference.XWorldLimits,...
%     fliplr(hdr.RasterReference.YWorldLimits));
% 
% cropr=floor(cropr);
% cropc=floor(cropc);

%allocate 3 & 4D cubes with measurements for all days
sz0=[hdr.RasterReference.RasterSize nbands length(matdates)];
refl=zeros(sz0);
SolarZenith=zeros([sz0(1) sz0(2) sz0(4)]);
cloudmask=false([sz0(1) sz0(2) sz0(4)]);
pxweights=zeros([sz0(1) sz0(2) sz0(4)]);

parfor i=1:length(matdates)
    isodate=datenum2iso(matdates(i),7);
    %allocate daily cubes
    refl_=NaN(sz);
    SolarZenith_=NaN([sz(1) sz(2)]);
    SolarAzimuth_=NaN([sz(1) sz(2)]);
    cloudmask_=false([sz(1) sz(2)]);
    pxweights_=zeros([sz(1) sz(2)]);
    %load up each tile
    for k=1:length(tiles)
        tile=tiles{k};
        
        %get topofile for each tile
        d=dir(fullfile(topodir,['*' tile '*.h5']));
        assert(~isempty(d),'%s empty\n',topodir);
        topofile=fullfile(topodir,d.name);
        
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
            
            [~,aWeights,~] = weightMOD09(f,topofile);
            pxweights_(r,c)= aWeights;
                        
%             for j=1:length(sunnames)
                x=single(GetMOD09GA(f,'SolarZenith'));
                if any(isnan(x(:)))
                    x = inpaint_nans(double(x),4);
                end
                x=imresize(x,tsiz(k,:));
                SolarZenith_(r,c)=imresize(x,tsiz(k,:));
                
%             end
                x=single(GetMOD09GA(f,'SolarAzimuth'));
                if any(isnan(x(:)))
                    x = inpaint_nans(double(x),4);
                end
                x=imresize(x,tsiz(k,:));
                SolarAzimuth_(r,c)=imresize(x,tsiz(k,:));

            %get all band reflectance
            sr=GetMOD09GA(f,'allbands');
            
            %create cloud mask
            S=GetMOD09GA(f,'state');
            mod09gacm=imresize(S.cloud==1,tsiz(k,:));
            swir=sr(:,:,swir_b);       
            cm = mod09gacm & swir > swir_cloud_thresh;
            
            refl_(r,c,:)=sr;
            cloudmask_(r,c)=cm;
            fprintf('loaded, corrected, and created masks for tile:%s date:%i \n',...
                tile,isodate);
        else
            fprintf('MOD09GA for %s %i not found,skipped\n',tile,isodate);
        end
    end
    
    %crop everything to topofile hdr
%     fn=fieldnames(o);
    
%     for j=1:length(fn)
%        if strcmp(fn{j},'cloudmask')
%            method='nearest';
%        else
%            method='linear';
%        end
       refl_=rasterReprojection(refl_,BigR,mstruct,...
            hdr.ProjectionStructure,'rasterref',hdr.RasterReference)
       cloudmask_=rasterReprojection(cloudmask_,BigR,mstruct,...
            hdr.ProjectionStructure,'rasterref',hdr.RasterReference,...
            'Method','nearest');
       SolarZenith_=rasterReprojection(SolarZenith_,BigR,mstruct,...
            hdr.ProjectionStructure,'rasterref',hdr.RasterReference);
       SolarAzimuth_=rasterReprojection(SolarAzimuth_,BigR,mstruct,...
            hdr.ProjectionStructure,'rasterref',hdr.RasterReference);
       pxweights_=rasterReprojection(pxweights_,BigR,mstruct,...
            hdr.ProjectionStructure,'rasterref',hdr.RasterReference);
%     end
    %correct reflectance
    refl_c_=normalizeReflectance(refl_,Slope,Aspect,SolarZenith_,...
        SolarAzimuth_);
    refl(:,:,:,i)=refl_c_;
    SolarZenith(:,:,i)=SolarZenith_;
    cloudmask(:,:,i)=cloudmask_;
    pxweights(:,:,i)=pxweights_;
end

sz=size(refl);
filledCube=NaN(sz);
mvec=reshape(mask,[sz(1)*sz(2) 1]);

for i=1:size(refl,3) % for each band
    tic;
    bandcube=squeeze(refl(:,:,i,:));
    bandcube(cloudmask)=NaN; %set all the clouds to NaN
    %fill datacube
    vec=reshape(bandcube,[sz(1)*sz(2) sz(4)])'; %col major (days x pixels)
    parfor j=1:size(vec,2)
        if ~mvec(j)
            v=vec(:,j)
            t=isnan(v);
            if nnz(~t) >= 2
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

    parfor j=1:size(XX,3)
        if any(isnan(XX(:,:,j)) & ~mask,'all') %fill any remaining NaNs 
            %that could  not be temporally interpolated spatially
             XX(:,:,j)=inpaint_nans(double(XX(:,:,j)),4);
        end
    end
    filledCube(:,:,i,:)=XX;
    t2=toc;
    fprintf('filled band:%i in %g min\n',i,t2/60);
end