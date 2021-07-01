function watermask=create_watermask(hdr,watermask_dir)
%get lat lon
        [lat,lon]=projinv(hdr.ProjectionStructure,...
            hdr.RasterReference.XWorldLimits,hdr.RasterReference.YWorldLimits);
        
        xtile=0:4:140;
        lontile=-180:10:170;
        
        ytile=0:4:52;
        lattile=80:-10:-60;
        
        y_idx=nan(4,1);
        x_idx=nan(4,1);
        
        latvec=round(lat(2)+5,-1):-10:round(lat(1)-5,-1);
        lonvec=round(lon(1)-5,-1):10:round(lon(2)+5,-1);
        
        for j=1:4
            if j<=length(latvec) && latvec(j) >= -60 && latvec(j) <=80
                y_idx(j)=find(latvec(j)==lattile);
            end
            if j<=length(lonvec) && lonvec(j) >= -180 && lonvec(j) <=170
                x_idx(j)=find(lonvec(j)==lontile);
            end
        end
        
        t= ~isnan(x_idx) & ~isnan(y_idx);
        y_idx=unique(y_idx(t));
        x_idx=unique(x_idx(t));
        
        tifnames=cell(length(y_idx),length(x_idx));
        for j=1:size(tifnames,1)
            for k=1:size(tifnames,2)
                tifnames{j,k}=fullfile(watermask_dir,...
                    sprintf('yearlyClassification2015-2020-%010.0f-%010.0f.tif',...
                    ytile(y_idx(j))*1e4,xtile(x_idx(k))*1e4));
            end
        end
        
        [Big,BigR]=mosaicTiles(tifnames(:),'tif','logical',int8(255));
        watermask=rasterReprojection(Big,BigR,...
            [],hdr.ProjectionStructure,'rasterref',hdr.RasterReference,...
            'Method','nearest');
end