function R0=createR0(loc,tile,matdates,sensor,net)
%create R0 (background snow-free reflectance) by examing a stack of images
%input :
%loc - directory where surface reflectance files live
%tile - e.g., 'h09v05' for MODIS
%matdates: dates to consider e.g. datenum([2019 8 1]):datenum([2019 9 30]);
%sensor: string w/ sensor name. Currently only works for 'MODIS'
%net - mccm seriesnetwork

switch sensor
    case 'MODIS'
        sz=[2400 2400 7];
        R0=nan([sz length(matdates)]);
        cloudOrsnowMask=false([sz(1) sz(2) length(matdates)]);
        SensorZenith=zeros(size(cloudOrsnowMask));
        parfor i=1:length(matdates)
            yr=year(matdates(i));
            doy=matdates(i)-datenum([yr 1 1])+1;
            searchName=sprintf("*%i%03i.%s*.hdf",yr,doy,tile);
            d=dir(fullfile(loc,searchName));
            if isempty(d)
                warning ('%s missing',searchName)
                continue
            elseif length(d) > 1
                error('multiple files found matching %s',searchName)
            else
                f=fullfile(d.folder,d.name);

                %get all band reflectance
                sr=GetMOD09GA(f,'allbands');
                sr(isnan(sr))=0;
                
                %create cloud mask
                z=GetMOD09GA(f,'SensorZenith');
                %expansive 
                SensorZenith(:,:,i)=imresize(z,[sz(1) sz(2)]);

                %create cloud mask
                S=GetMOD09GA(f,'state');
                %expansive 
                MOD35cm=imresize(S.cloud==1 | S.cloud==2,[sz(1) sz(2)]);
                
                %new MccM approach
                I=pxFeatures(sr,sz(3));
                %scale to integers for CNN
                scaleFactor=10000;
                I=int16(I.*scaleFactor);
                %GPU runs out of memory
                C = semanticseg(I,net,'ExecutionEnvironment','cpu');
                cm = C == 'cloud' ;

                %gets brightest part, but not full dry lake area
                saltpan=imresize(S.saltpan,[sz(1) sz(2)]);                
                
                cloudOrsnowMask(:,:,i) = cm | C=='snow' | MOD35cm & ~saltpan;
                
                R0(:,:,:,i)=sr;
                fprintf('loaded %s\n',f)
            end
        end
        
        %compute ndvi
        ndvi=squeeze((R0(:,:,2,:)-R0(:,:,1,:))./(R0(:,:,2,:)+R0(:,:,1,:)));
        %compute ndsi
        ndsi=squeeze((R0(:,:,4,:)-R0(:,:,6,:))./(R0(:,:,4,:)+R0(:,:,6,:)));
        %set clouds, snow to nan for ndvi only (not ndsi)
        ndvi(cloudOrsnowMask)=nan;
        
        %set off nadir shots to nan
        ndvi(SensorZenith>30)=nan;
        ndsi(SensorZenith>30)=nan;

        %compute max and max index for ndvi (max green up)
        [max_ndvi,idx_max_ndvi]=max(ndvi,[],3);

        %compute min and min index for ndsi (snow cover minimum)
        [min_ndsi,idx_min_ndsi]=min(ndsi,[],3);
        
        %compute b3 min and min ind
        b3=squeeze(R0(:,:,3,:));

        b3(b3<0.1)=nan; %too dark

        %compute min and min index for b3
        [min_b3,idx_min_b3]=min(b3,[],3);
                
        X=nan(sz);
        for i=1:sz(1)
            for j=1:sz(2)
                %spectra for max ndvi day
                max_ndvi_vals=squeeze(R0(i,j,:,idx_max_ndvi(i,j)));
                %min_ndsi_vals=squeeze(R0(i,j,:,idx_min_ndsi(i,j)));
                %spectra for min b3 day
                min_b3_vals=squeeze(R0(i,j,:,idx_min_b3(i,j)));
                
                if min_ndsi(i,j) > 0
                    %then use the spectra from min b3 day
                    X(i,j,:)=min_b3_vals;
                    %X(i,j,:)=min_ndsi_vals;
                else %use spectra from max ndvi day
                    X(i,j,:)=max_ndvi_vals;
                end
            end
        end

  R0=X;
  %interpolate NaNs
  for j=1:size(R0,3)
    R0(:,:,j)=inpaint_nans(double(squeeze(R0(:,:,j))),3);
  end
    otherwise
        error("only MODIS implemented")
end
end

