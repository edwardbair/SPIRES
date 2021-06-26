function out=run_spires_landsat(r0dir,rdir,demfile,Ffile,shade,tolval,...
    fsca_thresh,dust_rg_thresh,grain_thresh,dust_thresh,CCfile,...
    WaterMaskfile,CloudMaskfile,fIcefile,...
    el_cutoff,subset)

%run spires  for a landsat scene
% r0dir - R0 directory, must contain geotiff surface reflectances from USGS
% rdir - R directory
% demfile - matfile containing dem in m in same projection as R&R0
% must contain: Z, elevation in m;
% then hdr struct with fields: RefMatrix, RasterReference, and
% ProjectionStructure
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

solarZ=getOLIsolar(rdir);

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

solarZmat=ones(size(dem.Z)).*solarZ;

%get R0 refl and reproject to hdr
R0=getOLIsr(r0dir,dem.hdr);

nanmask=all(isnan(R0.bands),3);

%snow-covered scene and reproject to hdr
R=getOLIsr(rdir,dem.hdr);

%load adjustment files

adjust_vars={'cloudmask','fice','cc','watermask'};

for i=1:length(adjust_vars)
if i==1
    m=matfile(CloudMaskfile);  
elseif i==2
    m=matfile(fIcefile);
elseif i==3
    m=matfile(CCfile);
elseif i==4
    m=matfile(WaterMaskfile);
end

if ~isempty(regexp(adjust_vars{i},'.*mask.*','ONCE'))
    method='nearest';
    logicalflag=true;
else
    method='linear';
    logicalflag=false;
end
   
thdr=m.hdr;
% reproject if RefMatrices or sizes don't match
    if  any(dem.hdr.RefMatrix(:)~=thdr.RefMatrix(:)) || ...
                any(dem.hdr.RasterReference.RasterSize~=thdr.RasterReference.RasterSize)
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


m=nanmask | A.cloudmask | A.watermask;

o=run_spires(R0.bands,R.bands,solarZmat,Ffile,m,shade,...
    grain_thresh,dust_thresh,tolval,dem.hdr,red_b,swir_b,[]);

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
el_mask=dem.Z<el_cutoff;
ifsca(el_mask)=0;

% set pixels outside boundary, in cloudy mask, or in shade
ifsca(nanmask | A.cloudmask | A.watermask)=NaN;

igrainradius=single(o.grainradius);
igrainradius(isnan(ifsca) | ifsca==0)=NaN;

idust=single(o.dust);
idust(igrainradius<=dust_rg_thresh)=0;
idust(isnan(igrainradius))=NaN;
    
out.fsca_raw=fsca_raw;
out.fsca=ifsca;
out.fshade=fshade;
out.grainradius=igrainradius;
out.dust=idust;
out.watermask=A.watermask;
out.cloudmask=A.cloudmask;
out.nodatamask=nanmask;
out.cc=A.cc;
out.mu0=cosd(solarZ);

out.hdr=dem.hdr;

et=toc(t1);
fprintf('total elapsed time %4.2f min\n',et/60);
end