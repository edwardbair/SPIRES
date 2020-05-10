function [out,fname]=fill_and_run_modis(tiles,matdates,hdfbasedir,...
    topodir,topofile,mask,R0,Ffile,dust_thresh,dustmask,tolval,...
    outloc,nameprefix)

% fills input (mod09ga) and runs spires
%input:
% tiles - tilenames, e.g. 'h08v05'
% if tile is a cell vector, assumption is made to mosaic multiple tiles
% together then crop/reproject to topofile. watermask, R0,& dustmask will
% all be assumed to match topofile hdr for spatial info
% matdates - matdates for cube
% hdfbasedir - where the MOD09GA HDF files live
% must have sub directories that correspond to entries in tile, e.g. h08v04
% topodir - directory for h5 topo files from from consolidateTopography, 
% part of TopoHorizon that contain topofiles for each tile in tiles, 
% e.g. h09v05dem_463m_Topography.h5
% topofile- h5 file name from consolidateTopography, part of TopoHorizons
% watermask- logical mask w/ ones for pixels to exclude (like water)
% R0 - background image (MxNxb). Recommend using time-spaced smoothed
% cube from a month with minimum fsca and clouds, like August or September,
% then taking minimum of reflectance for each band (b)
% Ffile, location of griddedInterpolant object that produces 
% reflectances for each band
% with inputs: grain radius, dust, cosZ, i.e. the look up table, band
% dust_thresh: threshold value for dust/grain retrievls, e.g. 0.90
% dustmask: mask of loctations (MxN) where dust can be retireved (1) or not
% (0)
% tol val: threshold for uniquetol spectra, higher runs faster, 0 runs all pixesl
% scalar e.g. 0.05
% outloc: path to write output
% nameprefix - name prefix for outputs, e.g. Sierra

%output:
%   out:
%   fsca: MxNxd
%   grainradius: MxNxd
%   dust: MxNxd
%also writes 1 month .mat files with those outputs

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
    fillMODIScube(tiles,rundates,hdfbasedir,topodir,topofile,mask,swir_b);
    out=run_spires(R0,R,solarZ,Ffile,mask,dust_thresh,dustmask,tolval,...
        hdr,red_b,swir_b);
    fname=fullfile(outloc,[nameprefix datestr(rundates(1),'yyyymm') '.mat']);
    mfile=matfile(fname,'Writable',true);
    %fsca
    t=isnan(out.fsca);
    out.fsca_raw=uint8(out.fsca*100);
    out.fsca_raw(t)=intmax('uint8'); % 255 is NaN
    mfile.fsca_raw=out.fsca_raw;
    %grain radius
    t=t | out.fsca==0 | isnan(out.grainradius); %sometimes grain radius is NaN
    out.grainradius=uint16(out.grainradius);
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
t2=toc(t1);
fprintf('completed in %5.2f hr\n',t2/60/60);