t1=tic;
[R,refl,solarZ,cloudmask,snowmask,weights]=...
    smoothMODIScube(tile,matdates,hdfdir,topofile,snowR);
out=run_scagd_modis(R0,R,solarZ,F,watermask,fsca_thresh);
toc(t1);