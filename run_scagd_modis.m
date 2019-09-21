function out=run_scagd_modis(R0,R,solarZ,F,watermask,fsca_thresh,varargin)
% run LUT version of scagd for 4-D matrix R
% produces cube of: fsca, grain size (um), and dust concentration (by mass)
% no filters are applied (e.g. persistence)
% input:
% R0 - background image (MxNxb). Recommend using time-spaced smoothed
% cube from a month with minimum fsca and clouds, like August or September,
% then taking minimum of reflectance for each band (b)
% R - 4D cube of time-space smoothed reflectances (MxNxbxd) with
% dimensions: x,y,band,day
% solarZ: solar zenith angles (MxNxd) for R
% F: griddedInterpolant object that produces reflectances for each band
% with inputs: grain radius, dust, cosZ
% watermask: logical mask, true for water
% fsca_thresh: min fsca cutoff, scalar e.g. 0.15
% (optional) pshade: (photometric) shade spectra (bx1); reflectances
% corresponding to bands
%output:
%   out:
%   fsca: MxNxd
%   grainradius: MxNxd
%   dust: MxNxd

%use shade endmember
doShade=false;
pshade=[];

if nargin==7
    pshade=varargin{1};
    doShade=true;
end

sz=size(R);

fsca=zeros([sz(1)*sz(2) sz(4)]);
grainradius=NaN([sz(1)*sz(2) sz(4)]);
dust=NaN([sz(1)*sz(2) sz(4)]);

cosZ=cosd(solarZ);

cosZ=reshape(cosZ,[sz(1)*sz(2) sz(4)]);
R=reshape(R,[sz(1)*sz(2) sz(3) sz(4)]);
R0=reshape(R0,[sz(1)*sz(2) sz(3)]);
watermask=reshape(watermask,[sz(1)*sz(2) 1]);

veclen=sz(1)*sz(2);

for i=1:sz(4)
    thisR=squeeze(R(:,:,i));
    thiscosZ=squeeze(cosZ(:,i));
    tic;
    parfor j=1:veclen
        cZ=thiscosZ(j); %cosZ scalar (sometimes NaN on MOD09GA)
        wm=watermask(j); %watermask scalar
        if ~wm && ~isnan(cZ)
            pxR=squeeze(thisR(j,:)); %reflectance vector
            pxR0=squeeze(R0(j,:)); %background reflectance vector
            NDSI=(pxR(4)-pxR(6))/(pxR(4)+pxR(6));
            if NDSI > 0
                if doShade
                    o=speedyinvert(pxR,pxR0,cZ,F,pshade);
                    fsca(j,i)=o.x(1)/(1-o.x(2)); %normalize by fshade
                    grainradius(j,i)=o.x(3);
                    dust(j,i)=o.x(4);
                else
                    o=speedyinvert(pxR,pxR0,cZ,F);
                    fsca(j,i)=o.x(1);
                    grainradius(j,i)=o.x(2);
                    dust(j,i)=o.x(3);
                end
            end
        else
            fsca(j,i)=NaN;
        end
    end
    t2=toc;
    fprintf('done w/ day %i in %g min\n',i,t2/60);
end

fsca=reshape(fsca,[sz(1) sz(2) sz(4)]);

grainradius=reshape(grainradius,[sz(1) sz(2) sz(4)]);
dust=reshape(dust,[sz(1) sz(2) sz(4)]);

fsca(fsca<fsca_thresh)=0;
grainradius(fsca==0)=NaN;
dust(fsca==0)=NaN;

out.fsca=fsca;
out.grainradius=grainradius;
out.dust=dust;