tic;
[R,refl,solarZ,cloudmask,snowmask]=...
    smoothMODIScube(tile,matdates,hdfdir,topofile,snowR,fsca_nPersist);
tic;
out=run_scagd_modis(R0,R,solarZ,F,watermask,fsca_nPersist,...
    fsca_thresh,grainradius_nPersist,grainradius_thresh);
toc;