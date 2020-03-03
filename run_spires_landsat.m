function out=run_spires_landsat(r0dir,rdir,topofile,...
    Ffile,tolval,fsca_thresh,dust_thresh,pshade,CCfile,cloudmask,...
    subset)
%run spires  for a landsat scene
% r0date - date for background scene in yyyymmdd, e.g. 20180923
% r0dir - R0 directory, must contain geotiff surface reflectances from USGS
% rdate - date for snow scene in yyyymmdd, e.g. 20190131
% rdir - R directory
% topofile - h5 file containing scene topgraphy, outut of
% consolidate_topograpy
% Ffile - mat file containing F griddedInterpolant R=F(grain radius (um),...
% dust (ppmw),solarZenith (degrees),bands (scalar))
% tolval - uniquetol tolerance, e.g. 0.05 for separating unique spectra
% fsca_thresh - minumum fsca value for snow detection, values below are set to
% zero, e.g. 0.15, scalar
% dust_thresh - minumum fsca value for dust detection, pixels below are
% interpolated, e.g. 0.99, scalar
% pshade - physical shade endmember, vector, bandsx1
% CCfile - location of .mat
% cloudmask - cloudmask, logical
% canopy cover - canopy cover for pixels 0-1, size of scene
% also need RefMatrix and ProjectionStructure
% subset - either empty for none or [row1 row2;col1 col2], where are
% row1/col1 are the starting pixels and row2/col2 are the end pixels,
% e.g. for MMSA on p42r34, % [3280 3460;3740 3920]
%takes a while if not subsetting, e.g. p42r34 


%output
% o struct with fields:
% fsca,0-1, canopy cover adj
% grainradius, um
% dust,0-1
% shade,0-1
% all the size of the first two dimension of R or R0

red_b=3;
swir_b=6;

%do terrain first
%need to account for shaded pixels

%get sun data, R

t1=tic;

[solarZ,phi0]=getOLIsolar(rdir);

[Slope,hdr]=GetTopography(topofile,'slope');
Aspect=GetTopography(topofile,'aspect');

mu=sunslope(cosd(solarZ),180-phi0,Slope,Aspect);

h=GetHorizon(topofile,180-phi0);
%in sun if solarZ > 10 deg, shaded if solarZ <= 10 deg
smask= (90-acosd(mu))-h > 10;

%if crop w/o reprojection
if ~isempty(subset)
    rl=subset(1,1):subset(1,2);
    cl=subset(2,1):subset(2,2);
    [x,y]=pixcenters(hdr.RefMatrix,hdr.RasterReference.RasterSize,'makegrid');
    x=x(rl,cl);
    y=y(rl,cl);
    hdr.RefMatrix=makerefmat(x(1,1),y(1,1),x(1,2)-x(1,1),y(2,1)-y(1,1));
    hdr.RasterReference=refmatToMapRasterReference(hdr.RefMatrix,size(x));
    mu=mu(rl,cl);
    smask=smask(rl,cl);
    Slope=Slope(rl,cl);
    Aspect=Aspect(rl,cl);
    cloudmask=cloudmask(rl,cl);
end

smask=~imfill(~smask,8,'holes'); %fill in errant holes

%get R0 refl and reproject to hdr
R0=getOLIsr(r0dir,hdr);

nanmask=all(isnan(R0.bands),3);

%get sun data
[solarZR0,phi0R0]=getOLIsolar(r0dir);

%snow-covered scene and reproject to hdr
R=getOLIsr(rdir,hdr);

% load and reproject canopy data to hdr
CC=load(CCfile);
cc=rasterReprojection(CC.cc,CC.RefMatrix,CC.ProjectionStructure,...
    hdr.ProjectionStructure,'rasterref',hdr.RasterReference);
cc(isnan(cc))=0;
cc=double(cc);

t=normalizeReflectance(R.bands,Slope,Aspect,solarZ,phi0);
t0=normalizeReflectance(R0.bands,Slope,Aspect,solarZR0,phi0R0);

o=run_spires(t0,t,acosd(mu),Ffile,~smask | nanmask | cloudmask,...
    fsca_thresh,pshade,dust_thresh,tolval,cc,hdr,red_b,swir_b);

% spatial interpolation
ifsca=o.fsca;
ifsca(nanmask)=0;
ifsca(~smask & ~nanmask)=NaN;
ifsca=inpaint_nans(ifsca,4);

ifsca=ifsca./(1-cc);
ifsca(ifsca>1)=1;
ifsca(ifsca<fsca_thresh)=0;
ifsca(nanmask)=NaN;

igrainradius=o.grainradius;
igrainradius(nanmask)=0;
igrainradius(~smask & ~nanmask)=NaN;
igrainradius(igrainradius>1190)=NaN;
igrainradius=inpaint_nans(igrainradius,4);
igrainradius(ifsca==0)=NaN;
igrainradius(nanmask)=NaN;

idust=o.dust;
idust(nanmask)=0;
idust(~smask & ~nanmask)=NaN;
idust=inpaint_nans(idust,4);
idust(ifsca==0)=NaN;
idust(nanmask)=NaN;

out.fsca=ifsca;
out.grainradius=igrainradius;
out.dust=idust;
out.shade=o.shade;
out.hdr=hdr;

et=toc(t1);
fprintf('total elapsed time %4.2f min\n',et/60);

end