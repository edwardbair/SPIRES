function R0=createR0(loc,tile,matdates,sensor)
%create R0 (background snow-free reflectance) by examing a stack of images
%input :
%loc - directory where surface reflectance files live
%tile - e.g., 'h09v05' for MODIS
%matdates: dates to consider e.g. datenum([2019 8 1]):datenum([2019 9 30]);
%output R0
%sensor: string w/ sensor name. Currently only works for 'MODIS'

switch sensor
    case 'MODIS'
        R0=nan([2400 2400 7 length(matdates)]);
        for i=1:length(matdates)
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
                fname=fullfile(d.folder,d.name);
                R0(:,:,:,i)=GetMOD09GA(fname,'allbands');
                fprintf('loaded %s\n',fname)
            end
        end
        [~,idx]=min(squeeze(R0(:,:,3,:)),[],3);
        X=nan([size(idx) size(R0,3)]);
       sprintf('')
        for i=1:size(idx,1)
            for j=1:size(idx,2)
            X(i,j,:)=squeeze(R0(i,j,:,idx(i,j)));
            end
        end
  R0=X;
    otherwise
        error("only MODIS implemented")
end
end

