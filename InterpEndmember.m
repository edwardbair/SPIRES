function InterpEndmember(h5files,thresh_var,thresh,target_var,baseloc)
%spatially interpolate Endmember using IDW
%inputs:
%h5files: cell, Nx1 of h5 files to process,
%see GetEndmember for format
%thresh_var: variable for postings, e.g. 'snow'
%thresh: threshold value for postings, e.g. 0.95
%target_var: variable to interpolate, e.g. 'dust'
%baseloc: h5 location for read/write, e.g. '/Grid/MODIS_GRID_500m/'

for i=1:length(h5files)
    h5file=h5files{i};
    
    thisDataSet=[baseloc target_var];
    info = h5info(h5file,thisDataSet);

    divisor=h5readatt(h5file,thisDataSet,'divisor');
    switch info.Datatype.Type
        case 'H5T_STD_U8LE'
            dtype='uint8';
        case 'H5T_STD_U16LE'
            dtype='uint16';
    end
        
    [X,matdates,hdr]=GetEndmember(h5file,thresh_var);
    sz=size(X);
    Y=GetEndmember(h5file,target_var);
    [x,y]=pixcenters(hdr.RefMatrix,[sz(1) sz(2)],'makegrid');
    xv=x(:);
    yv=y(:);
    sz=size(X);
    for j=1:length(matdates)
        thisX=X(:,:,j);
        thisY=Y(:,:,j);
        thisX=thisX(:);
        thisY=thisY(:);
%         interpY=NaN(size(thisY));
        interpY=thisY(:);
        t=thisX>=thresh;
        if nnz(t) >= 3
            tt=~isnan(thisY) & thisX<thresh;
            idx=find(tt);
            iY=thisY(tt);
            tY=thisY(t);
            xvt=xv(t);
            yvt=yv(t);
            parfor k=1:length(idx)
                %inverse distance weighted average
                dist=pdist2([xvt yvt],[xv(idx(k)) yv(idx(k))]);
                w=1./dist;
                iY(k)=(sum(w.*tY))./(sum(w));
            end
            interpY(tt)=iY;
            interpY=reshape(interpY,[sz(1) sz(2)]);
            Y(:,:,j)=interpY;
        end
        fprintf('finished %s\n',datestr(matdates(j)));
    end
       
    Yint = float2integer(Y,divisor,0,dtype,0,max(Y(:)));
    h5write(h5file,[baseloc target_var],Yint)

    fprintf('wrote %s\n',h5file);
end



