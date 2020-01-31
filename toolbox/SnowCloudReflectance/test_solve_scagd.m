a=squeeze(R(r1,c1,:,1));
b=squeeze(R0(r1,c1,:,1));
c=squeeze(solarZ(r1,c1,1));

tic;o=speedyinvert(a,b,c,F,pshade,[]);toc;
tic;o=speedyinvert_mex(a,b,c,Fmat,Fkey,pshade,0,0);toc;

out=run_scagd_modis(R0,R,solarZ,Fmat,Fkey,watermask,fsca_thresh,pshade);