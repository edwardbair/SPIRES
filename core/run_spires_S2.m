function out=run_spires_S2(r0dir,rdir,demfile,Ffile,shade,tolval,...
    fsca_thresh,dust_rg_thresh,grain_thresh,dust_thresh,CCfile,...
    WaterMaskfile,CloudMaskfile,fIcefile,...
    el_cutoff,subset)

%run spires  for a Sentinel2 scene
% r0dir - R0 directory
% rdir - R directory
% demfile - matfile containing dem in m in same projection as R&R0
% must contain: Z, elevation in m;
% then RR, a mapcellref w/ ProjectedCRS field
% Ffile - mat file containing F griddedInterpolant R=F(grain radius (um),...
% dust (ppmw),solarZenith (degrees),bands (scalar))
% shade - shade endmeber, scalar or vector, length #bands
% tolval - uniquetol tolerance, e.g. 0.05 for separating unique spectra
% fsca_thresh - minumum fsca value for snow detection, values below are set to
% zero, e.g. 0.10, scalar
% dust_rg_thresh - min grain radius for dust calcl, e.g. 400 um
% grain_thresh - min fsca for grain size, e.g. 0.90;
% dust_thresh - min fsca for dust calc, e.g. 0.90;
% CCfile - location of .mat
% canopy cover - canopy cover for pixels 0-1, size of scene
% WaterMaskfile - water mask file location, contains watermask, logical
% 0 is no water, 1 is water
% CloudMaskfile - cloud mask file location, contains cloudmask, logical,
% 0 is no cloud, 1 is cloud
% fIcefile, fice file location, ice fraction 0-1, single, size of cloudmask
% also need RefMatrix and ProjectionStructure
% el_cutoff - elevation cutoff, m
% subset - either empty for none or [row1 row2;col1 col2], where are
% row1/col1 are the starting pixels and row2/col2 are the end pixels,

%output
% o struct with fields:
% fsca,0-1, canopy cover adj
% grainradius, um
% dust,0-1
% shade,0-1
% all the size of the first two dimension of R or R0

%note numbering corresponds to available bands & 
%how they are sorted by MATALB: 2,3,4,5,6,7,11,12,8a
red_b=3;
swir_b=9;
res='20'; %m resolution
  
t1=tic;

%get R0 refl 
R0=getS2sr(r0dir,res,subset);

%snow-covered scene
R=getS2sr(rdir,res,subset);

%load adjustment files
adjust_vars={'cloudmask','fice','cc','watermask','Z'};
%put ancillary data into struct A
for i=1:length(adjust_vars)
switch adjust_vars{i}
    case 'cloudmask'
    m=matfile(CloudMaskfile);  
    case 'fice'
    m=matfile(fIcefile);
    case 'cc'
    m=matfile(CCfile);
    case 'watermask'
    m=matfile(WaterMaskfile);
    case 'Z'
    m=matfile(demfile);
end
    A.(adjust_vars{i})=m.(adjust_vars{i});
end
A.RR=R.RR; %set mapcellsref from R 

%if subset on
if ~isempty(subset)
%if crop w/o reprojection
    for i=1:length(adjust_vars)
        X=A.(adjust_vars{i});
        X=X(subset(1,1):subset(1,2),subset(2,1):subset(2,2));
        A.(adjust_vars{i})=X;
    end
    A.RR=R0.RR; %set to cropped RR
end

%solar zenith
solarZmat=ones(size(A.Z)).*round(R.solarzenith);
%full mask
mask=A.cloudmask | A.watermask;

o=run_spires(R0.bands,R.bands,solarZmat,Ffile,mask,shade,...
    grain_thresh,dust_thresh,tolval,A.RR,red_b,swir_b,[]);

fsca_raw=single(o.fsca);
t0=fsca_raw==0; %track zeros to prevent 0/0 = NaN

%fshade correction
fshade=o.fshade;
A.cc(isnan(A.cc))=0;
A.fice(isnan(A.fice))=0;

ifsca=fsca_raw./(1-fshade-A.fice);
ifsca(ifsca>1 | ifsca<0)=1;
ifsca(t0)=0;

ifsca(ifsca<fsca_thresh)=0;

ifsca(ifsca>0 & A.cc>0)=1;

%elevation cutoff
el_mask=A.Z<el_cutoff;
ifsca(el_mask)=0;

% set masked values to nan
ifsca(A.cloudmask | A.watermask)=NaN;

igrainradius=single(o.grainradius);
igrainradius(isnan(ifsca) | ifsca==0)=NaN;

idust=single(o.dust);
idust(igrainradius<=dust_rg_thresh)=0;
idust(isnan(igrainradius))=NaN;
    
out.fsca_raw=fsca_raw;
out.fsca=ifsca;
out.fshade=single(fshade);
out.grainradius=igrainradius;
out.dust=idust;
out.watermask=A.watermask;
out.cloudmask=A.cloudmask;
out.nodatamask=isnan(fsca_raw);
out.cc=A.cc;
out.mu0=cosd(R.solarzenith);

out.RR=A.RR;

et=toc(t1);
fprintf('total elapsed time %4.2f min\n',et/60);
end