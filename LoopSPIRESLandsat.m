function LoopSPIRESLandsat(basedir,R0list,Rlist,Ffile,shade,tolval,...
    fsca_thresh,grain_thresh,dust_thresh,outdir,el_cutoff)
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
%"dustmask" in pPPPrRRR.mat with variable dustmask, logical
% R0list - directory list for R0, cell, Nx1
% Rlist - directory list for R, cell Nx1
% Ffile - mat file containing F griddedInterpolant R=F(grain radius (um),...
% dust (ppmw),solarZenith (degrees),bands (scalar))
% shade - shade endmeber, scalar or vector, length #bands
% tolval - uniquetol tolerance, e.g. 0.05 for separating unique spectra
% fsca_thresh - minumum fsca value for snow detection, values below are set to
% zero, e.g. 0.15, scalar
% grain_tresh - minimum fsca for grain size detection, e.g 0.25
% dust_thresh - minimum fsca value for dust & grain size detection e.g.
% 0.75
% outdir - where to write files out
% el_cutoff - elevation cutoff, m
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
    CloudMaskfile=fullfile(basedir,'wvmask','cloudmask',[Rlist{i},'_wvmask.mat']);
    cm=matfile(CloudMaskfile);
    hdr=cm.hdr;
    mask=cm.cloudmask;
    [x,y]=pixcenters(hdr.RefMatrix,hdr.RasterReference.RasterSize,'makegrid');
    x=x(~mask);
    y=y(~mask);
    [r,c]=map2pix(hdr.RefMatrix,x,y);
    subset=[min(r) max(r);min(c) max(c)];
    
    DustMaskfile=fullfile(basedir,'dustmask',fname);
    fIcefile=fullfile(basedir,'fice',fname);
   
    out=run_spires_landsat(r0dir,rdir,demfile,...
        Ffile,shade,tolval,fsca_thresh,grain_thresh,dust_thresh,...
        DustMaskfile,CCfile,...
        WaterMaskfile,CloudMaskfile,fIcefile,el_cutoff,subset);
    
    fn=fieldnames(out);
    outname=fullfile(outdir, [Rlist{i} '.mat']);
    m=matfile(outname,'Writable',true);
    
    for j=1:length(fn)
        m.(fn{j})=out.(fn{j});
    end
    
end