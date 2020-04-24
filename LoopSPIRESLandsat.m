function LoopSPIRESLandsat(basedir,R0list,Rlist,Ffile,tolval,...
    fsca_thresh,pshade,outdir,el_cutoff,subset)
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
% pshade - physical shade endmember, vector, bandsx1
% outdir - where to write files out
% el_cutoff - elevation cutoff, m
%subset - either empty for none or [row1 row2;col1 col2], where are
% row1/col1 are the starting pixels and row2/col2 are the end pixels,
% e.g. for MMSA on p42r34, % [3280 3460;3740 3920]
%note subset is based of DEM, as L8 has different sized scenes for
%different dates and everything is reprojected to match the dem
%takes a while if not subsetting, e.g. p42r34 

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
    WaterMaskfile=fullfile(basedir,'watermask',fname);
    CloudMaskfile=fullfile(basedir,'cloudmask',fname);
    DustMaskfile=fullfile(basedir,'dustmask',fname);
    fIcefile=fullfile(basedir,'fice',fname);
    
    out=run_spires_landsat(r0dir,rdir,demfile,...
        Ffile,tolval,fsca_thresh,DustMaskfile,pshade,CCfile,...
        WaterMaskfile,CloudMaskfile,fIcefile,el_cutoff,subset);
    
    fn=fieldnames(out);
    outname=fullfile(outdir, [Rlist{i} '.mat']);
    m=matfile(outname,'Writable',true);
    
    for j=1:length(fn)
        m.(fn{j})=out.(fn{j});
    end
    
end