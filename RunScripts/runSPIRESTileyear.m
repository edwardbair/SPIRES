function runSPIRESTileyear(code_dir,R0file,maskfile,topofile,...
    outloc,Ffile,hdfbasedir,tile,ccfile,ficefile,mccmfile,cores,WY,...
    b_R,dust_rg_thresh,...
    dust_thresh,el_cutoff,fsca_thresh,grain_thresh,maxdust,...
    maxgrainradius,mindust,mingrainradius,shade,tolval,windowSize,...
    windowThresh,Nd)
%for running each tile year from command prompt using batch mode

disp('loading files')
addpath(genpath(code_dir))

m=matfile(R0file);
R0=m.R0;

m=matfile(maskfile);
mask=m.mask;

m=matfile(mccmfile);
net=m.net;

m=matfile(ccfile);
cc=m.cc;

m=matfile(ficefile);
fice=m.fice;

if exist(outloc,'dir')==0
    mkdir(outloc);
end

disp('starting parpool')
parpool_check(cores)

matdates=datenum([WY-1 10 1]):datenum([WY 9 30]);

disp('running unmixing') 
fill_and_run_modis({tile},R0,matdates,...
hdfbasedir,topofile,mask,Ffile,shade,grain_thresh,dust_thresh,tolval,...
outloc,tile,net);

disp('smoothing')
smoothSPIREScube(tile,outloc,matdates,...
windowSize,windowThresh,mingrainradius,maxgrainradius,mindust,maxdust,...
mask,topofile,el_cutoff,fsca_thresh,cc,fice,b_R,dust_rg_thresh,true,Nd);
