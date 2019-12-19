function out=smooth_and_run(tile,matdates,hdfdir,topofile,watermask,...
    R0,F,pshade,fsca_thresh,outloc,grainradius_nPersist,el_cutoff)

% smooths input (mod09ga) and runs scagd, then smooths again
%input:
%tile - tilename, e.g. 'h08v05'
%matdates - matdates for cube
%hdfdir - where the MOD09GA HDF files live for a certain tile, e.g. h08v04
%topofile- h5 file name from consolidateTopography, part of TopoHorizons
%watermask- logical mask w/ ones for water
% R0 - background image (MxNxb). Recommend using time-spaced smoothed
% cube from a month with minimum fsca and clouds, like August or September,
% then taking minimum of reflectance for each band (b)

% F: griddedInterpolant object that produces reflectances for each band
% with inputs: grain radius, dust, cosZ, i.e. the look up table
% pshade: (photometric) shade spectra (bx1); reflectances
% corresponding to bands

% fsca_thresh: min fsca cutoff, scalar e.g. 0.15
% outloc: path to write output
% grainradius_nPersist: min # of consecutive days needed with normal 
%grain sizes to be kept as snow, e.g. 7
%el_cutoff, min elevation for snow, m - scalar, e.g. 1500
% dust_fsca_thresh - don't trust dust values with fsca below this, scalar,
% e.g. 0.85

%output:
%   out:
%   fsca: MxNxd
%   grainradius: MxNxd
%   dust: MxNxd
t1=tic;
%run in one month chunks and save output
%79 hr for 1 tile for 1 year
dv=datevec(matdates);
m=unique(dv(:,2),'stable');

for i=1:length(m)
    idx=dv(:,2)==m(i);
    rundates=matdates(idx);
    [R,~,solarZ,~,weights]=...
    smoothMODIScube(tile,rundates,hdfdir,topofile,watermask);
    out=run_scagd_modis(R0,R,solarZ,F,watermask,fsca_thresh,pshade);
    fname=fullfile(outloc,[datestr(rundates(1),'yyyymm') '.mat']);
    mfile=matfile(fname,'Writable',true);
    mfile.fsca=uint8(out.fsca*100);
    mfile.grainradius=uint16(out.grainradius);
    mfile.dust=uint16(out.dust*10);
    mfile.weights=uint8(weights*100);
    mfile.matdates=rundates;
    fprintf('wrote %s \n',fname);
end
%refilter and smooth
out=smoothSCAGDcube(outloc,matdates,...
    grainradius_nPersist,watermask,topofile,el_cutoff,fsca_thresh);

toc(t1)