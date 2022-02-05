function createTopoHorizonsMODIS(Z,ZR,tile,outloc)
%wrapper to create and save topohorizons for a tile 
%input: Z - global dem, e.g. GMTED2010
% ZR - rasterref w/ projCRS field for Z
% tile - tilename, char e.g. 'h24v05'
% outloc - directory to save output to. Note HD5 error for remote
% locations, need local scratch space

[R,p,RR] = sinusoidProjMODtile(tile);
ReferencingMatrix=R.RefMatrix_500m;
RR=RR.RasterReference_500m;
out.elevation=rasterReprojection(Z,ZR,'rasterref',RR);
out.elevation=single(out.elevation);
[A,H] = horizonAllDirections(out.elevation,RR);
out.viewfactor = viewFactor(A,H,out.elevation,RR);
out.viewfactor = round(out.viewfactor*100);
out.viewfactor = out.viewfactor / 100;
[out.slope,out.aspect] = topographicSlope(out.elevation,RR);
out.slope=round(out.slope);
out.aspect=floor(out.aspect);


%save
%horizons + proj information
fname=fullfile(outloc,sprintf('%sdem_463m_Topography.h5',tile));
if exist(fname,"file")
    delete(fname);
end
saveHorizon(fname,RR,A,H,'proj',p);

dnames={'elevation','slope','aspect','viewfactor'};
units={'m','degrees','degrees','degrees','unitless'};
dtypes={'int16','int16','int16','uint16'};
scalefactors=[1 1 1 100];
group='/Grid';
deflateLevel=9;

%still need RefMatrix for backward compatability
h5writeatt(fname,group,'ReferencingMatrix',ReferencingMatrix);

for i=1:length(dnames)
     dname=dnames{i};
%     switch dname
%         case 'elevation'
%             scaleFactor = 1;
%         case 'viewfactor'
%             scaleFactor = single(intmax(dtypes{i}));
%         otherwise
%             scaleFactor=single(intmax(dtypes{i}))/(max(out.(dname)(:)-min(out.(dname)(:))));
%     end
%     o = float2integer(out.(dname),scaleFactor,0,dtypes{i});
    t=isnan(out.(dname));
    o=cast(double(out.(dname)*scalefactors(i)),dtypes{i});
    dfill = intmin(dtypes{i});
    o(t)=dfill;

    h5create(fname,[group '/' dname],size(o),...
        'Deflate',deflateLevel,...
        'ChunkSize',[size(o,1) size(o,2)],...
        'FillValue',dfill,...
        'DataType',dtypes{i});
    h5write(fname,[ group '/' dname],o)
    h5writeatt(fname,[ group '/' dname],'divisor',scalefactors(i));
    h5writeatt(fname,[ group '/' dname],'units',units{i});
end

