function out=run_scagd(R0,R,solarZ,F,watermask)
%run LUT version of scagd for 4-D matrix R
%produces cube of: fsca, grain size (um), and dust concentration (by mass)
%input:
% R0 - background image (MxNxb). Recommend using time-spaced smoothed 
% cube from a month with minimum fsca and clouds, like August or September, 
% then taking minimum of reflectance for each band (b)
% R - 4D cube of time-space smoothed reflectances (MxNxbxd) with
% dimensions: x,y,band,day
% solarZ: solar zenith angles (MxNxd) for R
% F: griddedInterpolant object that produces reflectances for each band 
% with inputs: grain radius, dust, cosZ
% pshade: photometric shade spectra (bx1); reflectances 
% corresponding to bands
%output: 
%   out:
%   fsca: MxNxd
%   grainradius: MxNxd
%   dust: MxNxd

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
    parfor j=1:veclen
        if ~watermask(j)
            pxR=squeeze(thisR(j,:));
            pxR0=squeeze(R0(j,:));
            cZ=squeeze(thiscosZ(j));
            NDSI=(pxR(4)-pxR(6))/(pxR(4)+pxR(6));
                if NDSI > 0 && ~isnan(cZ)
                    % make sure all inputs are scalars or 
                    o=speedyinvert(pxR,pxR0,cZ,F);
                    %fsca_n=o.x(1);
%                     %if fsca+fshade > 1
%                     if o.x(1)+o.x(2) > 1
%                         fsca_n=1; %assume fsca=1
%                     else
%                         fsca_n=o.x(1)/(1-o.x(2)); %normalize by fshade
%                     end
                    %fsca_n(fsca_n>1)=1;
                    %fsca(j,i)=fsca_n;
                    fsca(j,i)=o.x(1);
                    grainradius(j,i)=o.x(2);
%                     grainradius(j,i)=o.x(2);
                    dust(j,i)=o.x(3);
%                     dust(j,i)=o.x(3);
                end
        else
            fsca(j,i)=NaN;
        end
    end
end


fsca=reshape(fsca,[sz(1) sz(2) sz(4)]);
grainradius=reshape(grainradius,[sz(1) sz(2) sz(4)]);
dust=reshape(dust,[sz(1) sz(2) sz(4)]);

nPersist=5;
thresh=0.15;
    
fsca=snowPersistenceFilter(fsca,nPersist,thresh);
fsca(fsca<thresh)=0;
grainradius(fsca==0)=NaN;
dust(fsca==0)=NaN;

out.fsca=fsca;
out.grainradius=grainradius;
out.dust=dust;

