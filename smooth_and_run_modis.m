function out=smooth_and_run_modis(tile,matdates,hdfdir,topofile,watermask,...
    R0,Ffile,pshade,dust_thresh,tolval,fsca_thresh,outloc,...
    grainradius_nPersist,el_cutoff,cc,fice)

% smooths input (mod09ga) and runs spires, then smooths again
%input:
%tile - tilename, e.g. 'h08v05'
%matdates - matdates for cube
%hdfdir - where the MOD09GA HDF files live for a certain tile, e.g. h08v04
%topofile- h5 file name from consolidateTopography, part of TopoHorizons
%watermask- logical mask w/ ones for water
% R0 - background image (MxNxb). Recommend using time-spaced smoothed
% cube from a month with minimum fsca and clouds, like August or September,
% then taking minimum of reflectance for each band (b)

% Ffile, location of griddedInterpolant object that produces 
% reflectances for each band
% with inputs: grain radius, dust, cosZ, i.e. the look up table, band
% pshade: (photometric) shade spectra (bx1); reflectances
% corresponding to bands
% dust thresh: min value for dust retrievals, scalar e.g. 0.85
% tol val: threshold for uniquetol spectra, higher runs faster, 0 runs all pixesl
% scalar e.g. 0.05
% fsca_thresh: min fsca cutoff, scalar e.g. 0.15
% outloc: path to write output
% grainradius_nPersist: min # of consecutive days needed with normal 
% grain sizes to be kept as snow, e.g. 7
% el_cutoff, min elevation for snow, m - scalar, e.g. 1500
% cc - static canopy cover, single or double, same size as watermask,
% 0-1 for viewable gap fraction correction
% fice - ice fraction, single or double 0-1
%5.77 hr for h08v05 WY 2017 w/ 62 cores
%output:
%   out:
%   fsca: MxNxd
%   grainradius: MxNxd
%   dust: MxNxd

red_b=3;
swir_b=6;

t1=tic;
%run in one month chunks and save output
dv=datevec(matdates);
m=unique(dv(:,2),'stable');
[~,hdr]=GetTopography(topofile,'elevation');

for i=1:length(m)
    idx=dv(:,2)==m(i);
    rundates=matdates(idx);
    [R,~,solarZ,~,weights]=...
    smoothMODIScube(tile,rundates,hdfdir,topofile,watermask);
    out=run_spires(R0,R,solarZ,Ffile,watermask,fsca_thresh,pshade,...
        dust_thresh,tolval,cc,hdr,red_b,swir_b);
    fname=fullfile(outloc,[tile datestr(rundates(1),'yyyymm') '.mat']);
    mfile=matfile(fname,'Writable',true);
    %fsca
    t=isnan(out.fsca);
    out.fsca=uint8(out.fsca*100);
    out.fsca(t)=intmax('uint8'); % 255 is NaN
    mfile.fsca=out.fsca;
    %grain radius
    out.grainradius=uint16(out.grainradius);
    t=t | out.fsca==0;
    out.grainradius(t)=intmax('uint16'); %65535
    mfile.grainradius=out.grainradius;
    %dust
    out.dust=uint16(out.dust*10);
    out.dust(t)=intmax('uint16'); %65535
    mfile.dust=out.dust;
    %weights
    t=isnan(weights);
    out.weights=uint8(weights*100);
    out.weights(t)=intmax('uint8'); %255
    mfile.weights=out.weights;
    
    mfile.matdates=rundates;
    fprintf('wrote %s \n',fname);
end
%refilter and smooth
out=smoothSCAGDcube(tile,outloc,matdates,...
    grainradius_nPersist,watermask,topofile,el_cutoff,fsca_thresh,cc,fice);

%write out h5 cubes
fname=fullfile(outloc,[tile datestr(matdates(end),'yyyy') '.h5']);
if exist(fname,'file')
    delete(fname); 
end
fn={'fsca','grainradius','dust'};
fntarget={'snow_fraction','grain_size','dust'};
dtype={'uint8','uint16','uint16'};
divisors=[100 1 10];

for i=1:length(fn)   
    member=fntarget{i};
    Value=out.(fn{i});
    dS.(member).divisor=divisors(i);
    dS.(member).dataType=dtype{i};
    dS.(member).maxVal=max(Value(:));
    dS.(member).FillValue=intmax(dS.(member).dataType);
    writeh5stcubes(fname,dS,out.hdr,out.matdates,member,Value);
end
t2=toc(t1);
fprintf('completed in %5.2f hr\n',t2/60/60);