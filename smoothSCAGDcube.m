function out=smoothSCAGDcube(outloc,matdates,...
    grainradius_nPersist,watermask,topofile,el_cutoff,fsca_thresh)

for i=1:length(matdates)
    dv=datevec(matdates(i));
    fname=fullfile(outloc,[datestr(dv,'yyyymm') '.mat']);
    m=matfile(fname);
    if i==1
        fsca=zeros([size(m.fsca,1) size(m.fsca,2) length(matdates)],'single');
        weights=zeros([size(m.fsca,1) size(m.fsca,2) length(matdates)],'single');
        %note zeros due to uint storage
        grainradius=zeros([size(m.fsca,1) size(m.fsca,2) length(matdates)],'single');
        dust=zeros([size(m.fsca,1) size(m.fsca,2) length(matdates)],'single');    
    end
    ind=datenum(dv)-datenum([dv(1:2) 1])+1;%1st day of month
    
    %take careful note of scaling coefficients
    %fsca
    fsca_t=single(m.fsca(:,:,ind));
    t=fsca_t==intmax('uint8');
    fsca_t(t)=NaN;
    fsca_t=fsca_t/100;
    %weights    
    weights_t=single(m.weights(:,:,ind));
    t=weights_t==intmax('uint8');
    weights_t(t)=NaN;
    weights_t=weights_t/100;
    %grain radius
    grainradius_t=single(m.grainradius(:,:,ind));
    t=grainradius_t==intmax('uint16');
    grainradius_t(t)=NaN;
    %dust
    dust_t=single(m.dust(:,:,ind));
    t=dust_t==intmax('uint16');
    dust_t(t)=NaN;
    dust_t=dust_t/10;
    
    fsca(:,:,i)=fsca_t;
    weights(:,:,i)=weights_t;
    grainradius(:,:,i)=grainradius_t;
    dust(:,:,i)=dust_t;
end

Z=GetTopography(topofile,'elevation');
Zmask=Z < el_cutoff;
Zmask=repmat(Zmask,[1 1 length(matdates)]);

wm=repmat(watermask,[1 1 size(fsca,3)]);

fsca(Zmask | wm) = 0;

%create mask for cube where radius is > 50 & radius < 1190 for 7 or more days

gmask=snowPersistenceFilter(grainradius > 50 & grainradius < 1190,...
    grainradius_nPersist,1);

% set to NAN days that aren't in that mask but are not zero fsca
fsca(~gmask & ~fsca==0)=NaN;

newweights=weights;
newweights(isnan(fsca))=0;
%fill in and smooth NaNs
fsca=smoothDataCube(fsca,newweights,'mask',~watermask);
%get some small fsca values from smoothing - set to zero
fsca(fsca<fsca_thresh)=0;
fsca(wm)=NaN;

%create mask of any fsca for interpolation
anyfsca=any(fsca,3);

%convert zeros(uint) to back to NaN
%grainradius(grainradius==0)=NaN;
% create mask for low fsca (includes zeros)
lowfscamask= fsca > 0 & fsca < 0.30;
%set all low fsca values to NaN
grainradius(lowfscamask)=NaN;
% set all weights for NaNs to 0
newweights=weights;
newweights(isnan(grainradius))=0;

grainradius=smoothDataCube(grainradius,newweights,'mask',anyfsca);
grainradius(fsca==0 | isnan(fsca))=NaN;

%convert zeros back to NaN
%dust(dust==0)=NaN;
%set all dust values to NaN
%dust(lowfscamask)=NaN;
%set weights
%newweights=weights;
%newweights(isnan(dust))=0; %if its NaN its zero
%dust(isnan(dust))=0; % need zeros for spline smoothing 

dust=smoothDataCube(dust,newweights,'mask',anyfsca);
%double smoothing needed to get rid of bloc proc artifacts
%dust=smoothDataCube(dust,newweights,...
 %   'mask',anyfsca,'method','smoothingspline');
dust(fsca==0 | isnan(fsca))=NaN;


[~,hdr]=GetTopography(topofile,'elevation');

out.matdates=matdates;
out.hdr=hdr;
out.fsca=fsca;
out.grainradius=grainradius;
out.dust=dust;

end

