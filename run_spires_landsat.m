function out=run_spires_landsat(r0dir,rdir,demfile,Ffile,tolval,...
    fsca_thresh,DustMaskfile,CCfile,WaterMaskfile,CloudMaskfile,...
    fIcefile,el_cutoff,subset)

%run spires  for a landsat scene
% r0date - date for background scene in yyyymmdd, e.g. 20180923
% r0dir - R0 directory, must contain geotiff surface reflectances from USGS
% rdate - date for snow scene in yyyymmdd, e.g. 20190131
% rdir - R directory
% demfile - matfile containing dem in m in same projection as R&R0
% must contain: Z, elevation in m;
% then hdr struct with fields: RefMatrix, RasterReference, and
% ProjectionStructure
% Ffile - mat file containing F griddedInterpolant R=F(grain radius (um),...
% dust (ppmw),solarZenith (degrees),bands (scalar))
% tolval - uniquetol tolerance, e.g. 0.05 for separating unique spectra
% fsca_thresh - minumum fsca value for snow detection, values below are set to
% zero, e.g. 0.15, scalar
% DustMaskfile - dust mask file location, locations where dust can be
% estimated
% watermask
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
% e.g. for MMSA on p42r34, % [3280 3460;3740 3920]
%note subset is based of DEM, as L8 has different sized scenes for
%different dates and everything is reprojected to match the dem
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

t1=tic;

[solarZ,phi0]=getOLIsolar(rdir);

dem=load(demfile);

if ~isempty(subset)
%if crop w/o reprojection
    rl=subset(1,1):subset(1,2);
    cl=subset(2,1):subset(2,2);
    [x,y]=pixcenters(dem.hdr.RefMatrix,dem.hdr.RasterReference.RasterSize,...
        'makegrid');
    x=x(rl,cl);
    y=y(rl,cl);
    dem.hdr.RefMatrix=makerefmat(x(1,1),y(1,1),x(1,2)-x(1,1),y(2,1)-y(1,1));
    dem.hdr.RasterReference=refmatToMapRasterReference(dem.hdr.RefMatrix,...
        size(x));
    dem.Z=dem.Z(rl,cl);
end
    tt1=tic;
    [Slope,Aspect] = SlopeAzmProjected(dem.hdr.RefMatrix, ...
    dem.hdr.ProjectionStructure, dem.Z);
    [x,y]=pixcenters(dem.hdr.RefMatrix,size(Slope),'makegrid');
    [lat,lon]=minvtran(dem.hdr.ProjectionStructure,x,y);
    
    sinF=Horizons2Directions('earth',180-phi0,lat,lon,dem.Z);
    tt2=toc(tt1);
    fprintf('Horizon calc in %4.2f min\n',tt2/60);

    h=asind(sinF);
    mu=sunslope(cosd(solarZ),180-phi0,Slope,Aspect);

    %in sun if solarZ > 10 deg, shaded if solarZ <= 10 deg
    smask= (90-acosd(mu))-h > 10;

smask=~imfill(~smask,8,'holes'); %fill in errant holes

%get R0 refl and reproject to hdr
R0=getOLIsr(r0dir,dem.hdr);

nanmask=all(isnan(R0.bands),3);

%get sun data
[solarZR0,phi0R0]=getOLIsolar(r0dir);

%snow-covered scene and reproject to hdr
R=getOLIsr(rdir,dem.hdr);

%load adjustment files

adjust_vars={'cloudmask','fice','cc','watermask','dustmask'};

for i=1:length(adjust_vars)
if i==1
    m=matfile(CloudMaskfile);  
elseif i==2
    m=matfile(fIcefile);
elseif i==3
    m=matfile(CCfile);
elseif i==4
    m=matfile(WaterMaskfile);
elseif i==5
    m=matfile(DustMaskfile);
end

if ~isempty(regexp(adjust_vars{i},'.*mask.*','ONCE'))
    method='nearest';
    logicalflag=true;
else
    method='linear';
    logicalflag=false;
end
   
thdr=m.hdr;
% reproject if RefMatrices don't match
    if  any(dem.hdr.RefMatrix(:)~=thdr.RefMatrix(:))
        A.(adjust_vars{i})=rasterReprojection(m.(adjust_vars{i}),...
        thdr.RefMatrix,thdr.ProjectionStructure,...
        dem.hdr.ProjectionStructure,'rasterref',...
        dem.hdr.RasterReference,'Method',method);
    else
        A.(adjust_vars{i})=m.(adjust_vars{i});
    end
    if ~logicalflag
        A.(adjust_vars{i})=double(A.(adjust_vars{i}));
        t=isnan(A.(adjust_vars{i}));
        A.(adjust_vars{i})(t)=0;
    else
        A.(adjust_vars{i})=logical(A.(adjust_vars{i}));
    end
end

%normalizeReflectance
t=normalizeReflectance(R.bands,Slope,Aspect,solarZ,phi0);
t0=normalizeReflectance(R0.bands,Slope,Aspect,solarZR0,phi0R0);

m=~smask | nanmask | A.cloudmask | A.watermask;

o=run_spires(t0,t,acosd(mu),Ffile,m,A.dustmask,tolval,...
    dem.hdr,red_b,swir_b);

fsca_raw=single(o.fsca);
t0=fsca_raw==0 | A.cc==1; %track zeros to prevent 0/0 = NaN

%viewable gap correction
A.cc(isnan(A.cc))=0;
ifsca=fsca_raw./(1-A.cc);

% fice correction
A.fice(isnan(A.fice))=0;
ifsca=ifsca./(1-A.fice);

%have to copy matrix across 3rd dim to set to min(fsca,ice)
fice_r=repmat(A.fice,[1 1 size(ifsca,3)]);
t=ifsca<fice_r;
ifsca(t)=fice_r(t);

ifsca(ifsca>1)=1;
ifsca(t0)=0;
ifsca(ifsca<fsca_thresh)=0;

%elevation cutoff
el_mask=dem.Z<el_cutoff;
ifsca(el_mask)=0;

% set pixels outside boundary, in cloudy mask, or in shade
ifsca(nanmask | A.cloudmask | ~smask | A.watermask)=NaN;

igrainradius=single(o.grainradius);
igrainradius(isnan(ifsca) | ifsca==0)=NaN;

idust=single(o.dust);
idust(isnan(ifsca) | ifsca==0)=NaN;

out.fsca_raw=fsca_raw;
out.fsca=ifsca;
out.grainradius=igrainradius;
out.dust=idust;
out.watermask=A.watermask;
out.shademask=~smask;
out.cloudmask=A.cloudmask;
out.nodatamask=nanmask;
out.cc=A.cc;

out.hdr=dem.hdr;

et=toc(t1);
fprintf('total elapsed time %4.2f min\n',et/60);
end