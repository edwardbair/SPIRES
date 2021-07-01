function [Z,status]=create_srtm(hdr,basedir,srtmdir)
status=false;
        [lat,lon]=projinv(hdr.ProjectionStructure,...
            hdr.RasterReference.XWorldLimits,hdr.RasterReference.YWorldLimits);
        latvals=floor(lat(1)):floor(lat(2));
        lonvals=floor(lon(1)):floor(lon(2));
        hgtnames=cell(length(latvals),length(lonvals));
        
        for j=1:length(latvals)
            for k=1:length(lonvals)
                if latvals(j) < 0
                    hemi='S';
                else
                    hemi='N';
                end
                if lonvals(k) < 0
                    mer='W';
                else
                    mer='E';
                end
                srtmfname=sprintf('%s%02.0f%s%03.0f.SRTMGL1.hgt.zip',...
                    hemi,abs(latvals(j)),mer,abs(lonvals(k)));
                srtmname=fullfile(srtmdir,srtmfname);
                if exist(srtmname,'file') == 2
                    %use tempname to prevent write collisions
                    tmpName = tempname;
                    tmpdir=fullfile(basedir,'temp',tmpName);
                    o=unzip(srtmname,tmpdir);
                    hgtnames{j,k}=o{1};
                else
                    fprintf('couldnt find %s \n',srtmname);
                end
            end
        end
        hgtnames=hgtnames(:);
        
        %remove empty names
        t=cellfun(@isempty,hgtnames);
        hgtnames(t)=[];
        try
            [Big,BigR]=mosaicTiles(hgtnames,'hgt','int16',intmin('int16'));
            Z=rasterReprojection(single(Big),BigR,[],hdr.ProjectionStructure,...
                'rasterref',hdr.RasterReference);
        catch
            fprintf('error w/ mosaicing %s\n',R_ID);
            Z=NaN(size(sr.bands,1),size(sr.bands,2));
            status=true;
        end
end