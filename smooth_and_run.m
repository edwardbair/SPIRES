function smooth_and_run(tile,matdates,hdfdir,topofile,watermask,...
    R0,F,fsca_thresh,outloc)
t1=tic;
%run in one month chunks and save output
dv=datevec(matdates);
m=unique(dv(:,2),'stable');

for i=1:length(m)
    idx=dv(:,2)==m(i);
    rundates=matdates(idx);
    [R,~,solarZ,~,weights]=...
    smoothMODIScube(tile,rundates,hdfdir,topofile,watermask);
    out=run_scagd_modis(R0,R,solarZ,F,watermask,fsca_thresh);
    fname=fullfile(outloc,[datestr(rundates(1),'yyyymm') '.mat']);
    mfile=matfile(fname,'Writable',true);
    mfile.fsca=uint8(out.fsca*100);
    mfile.grainradius=uint16(out.grainradius);
    mfile.dust=single(out.dust*10^11);
    mfile.weights=uint8(weights*100);
    mfile.matdates=rundates;
    fprintf('wrote %s \n',fname);
end

toc(t1)
%then pass output to a function that concatenates annual cubes,
%puts NaNs in for suspicious answers near limits, applies temporal
%filter, then smoothes concatenated cube