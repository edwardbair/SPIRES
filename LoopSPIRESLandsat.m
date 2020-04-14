function LoopSPIRESLandsat(basedir,R0list,Rlist,Ffile,tolval,...
    fsca_thresh,dust_thresh,pshade,outdir)
% call SPIRES Landsat in a loop
%input:
%basedir - base dir where L8 inputs live,
%must have subdirs:
%"cc" - canopy cover in pPPPrRRR.mat, where PPP and RRR are 3 digit path &
%rows, each file must contain cc as a canopy fraction, single
%"dem" - dem in pPPPrRRR_dem.mat with variable Z, single
%"fice" - ice fraction in pPPPrRRR.mat with variable fice, single
%"watermask" - watermask in pPPPrRRR.mat with variable watermask, logical
%"cloudmask" - cloudmask in pPPPrRRR.mat with variable cloudmask, logical
%"sr" - surface reflectance dir w/ subdirs in R & R0 list

% R0list - directory list for R0, cell, Nx1
% Rlist - directory list for R, cell Nx1
% Ffile - mat file containing F griddedInterpolant R=F(grain radius (um),...
% dust (ppmw),solarZenith (degrees),bands (scalar))
% tolval - uniquetol tolerance, e.g. 0.05 for separating unique spectra
% fsca_thresh - minumum fsca value for snow detection, values below are set to
% zero, e.g. 0.15, scalar
% dust_thresh - minumum fsca value for dust detection, pixels below are
% interpolated, e.g. 0.99, scalar
% pshade - physical shade endmember, vector, bandsx1
% outdir - where to write files out

for i=1:length(Rlist)
    rdir=fullfile(basedir,'sr',Rlist{i});
    r0dir=fullfile(basedir,'sr',R0list{i});
    [~,fpart]=fileparts(rdir);
    s=strsplit(fpart,'_');
    pathrow=['p' s{3}];
    pathrow=insertAfter(pathrow,4,'r');
    demfile=fullfile(basedir,'dem',[pathrow,'_dem.mat']);
    fname=[pathrow,'.mat'];
    CCfile=fullfile(basedir,'cc',fname);
    CloudMaskfile=fullfile(basedir,'cloudmask',fname);
    fIcefile=fullfile(basedir,'fice',fname);
    subset=[500 600; 500 600];
    
    out=run_spires_landsat(r0dir,rdir,demfile,...
        Ffile,tolval,fsca_thresh,dust_thresh,pshade,CCfile,...
        CloudMaskfile,fIcefile,subset);
    
    fn=fieldnames(out);
    outname=fullfile(outdir,fname);
    m=matfile(outname,'Writable',true);
    
    for j=1:length(fn)
        m.(fn{j})=out.(fn{j});
    end
    
end