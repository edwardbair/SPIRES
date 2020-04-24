function out=smoothSPIREScube(tile,outloc,matdates,...
    grainradius_nPersist,watermask,topofile,el_cutoff,fsca_thresh,cc,fice,...
    endconditions)
%function to smooth cube after running through SPIRES
% tile - tile name, string e.g. 'h08v05'
% outloc - output location, string
% matdates - datenum vector for image days
% grainradius_nPersist: min # of consecutive days needed with normal 
% grain sizes to be kept as snow, e.g. 7
% watermask- logical mask w/ ones for water 
% topofile- h5 file name from consolidateTopography, part of TopoHorizons
% el_cutoff, min elevation for snow, m - scalar, e.g. 1500
% fsca_thresh: min fsca cutoff, scalar e.g. 0.15
% cc - static canopy cover, single or doube, same size as watermask,
% 0-1 for viewable gap fraction correction
% fice - fraction of ice/neve, single or double, 0-1, mxn
% endconditions - string, end condition for splines for dust and grain size, 
% e.g. 'estimate' or 'periodic', see slmset.m

%output: struct out w/ fields
%fsca, grainradius, dust, and hdr (geographic info)

%1.34 hr/yr for h08v05, 2017
fprintf('reading %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));
for i=1:length(matdates)
    dv=datevec(matdates(i));
    fname=fullfile(outloc,[tile datestr(dv,'yyyymm') '.mat']);
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

fprintf('finished reading %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

[Z,hdr]=GetTopography(topofile,'elevation');
Zmask=Z < el_cutoff;
Zmask=repmat(Zmask,[1 1 length(matdates)]);

wm=repmat(watermask,[1 1 size(fsca,3)]);
fice=repmat(fice,[1 1 size(fsca,3)]);
fsca(Zmask | wm) = 0;

%create mask for cube where radius is > 50 & radius < 1190 for 7 or more days

gmask=snowPersistenceFilter(grainradius > 50 & grainradius < 1190,...
    grainradius_nPersist,1);

% set to NAN days that aren't in that mask but are not zero fsca
fsca(~gmask & ~(fsca==0))=NaN;

newweights=weights;
newweights(isnan(fsca))=0;
%fill in and smooth NaNs

fprintf('smoothing fsca %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

fsca=smoothDataCube(fsca,newweights,'mask',~watermask,...
   'method','smoothingspline','SmoothingParam',0.02);
%get some small fsca values from smoothing - set to zero
fsca(fsca<fsca_thresh)=0;
fsca(wm)=NaN;
%viewable gap correction
cc(isnan(cc))=0;
fsca=fsca./(1-cc);
fsca(fsca>1)=1;
% fice correction
fsca=fsca./(1-fice);
fsca(fsca>1)=1;
%min value of fsca is fice
t=fsca<fice;
fsca(t)=fice(t);

fprintf('finished smoothing fsca %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

fprintf('smoothing grain radius %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

%create mask of any fsca for interpolation
anyfsca=any(fsca,3);

% set all weights for NaNs to 0
newweights=weights;
newweights(isnan(fsca) | fsca==0)=0;

% grainradius=smoothDataCube(grainradius,newweights,'mask',anyfsca,...
%      'method','smoothingspline','SmoothingParam',0.75);
dF=cat(3,zeros(size(fsca,1,2)),diff(fsca,1,3));
%send logical cube for decreasing fsca
fcube=dF<=0;
grainradius=smoothDataCube(grainradius,newweights,'mask',anyfsca,...
    'method','slm','monotonic','increasing','fcube',fcube,'knots',-3,...
    'endconditions',endconditions);


grainradius(fsca==0 | isnan(fsca))=NaN;

fprintf('finished smoothing grain radius %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

%use a different approach
%replace high values with NaN
% parfor i=1:size(dust,3)
%    xx=squeeze(dust(:,:,i));
%    v=quantile(xx(:),0.95);
%    xx(xx>v)=NaN;
%    dust(:,:,i)=xx;
% end

fprintf('smoothing dust %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

%compute grain sizes differences
dG=cat(3,zeros(size(grainradius,1,2)),diff(grainradius,1,3));
%send logical cube for increasing grain sizes
fcube=dG>=0;

dust=smoothDataCube(dust,newweights,'mask',anyfsca,...
    'method','slm','monotonic','increasing','fcube',fcube,'knots',-3,...
    'endconditions',endconditions);
dust(fsca==0 | isnan(fsca))=NaN;

fprintf('finished smoothing dust %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

fprintf('writing cube %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

out.matdates=matdates;
out.hdr=hdr;
out.fsca=fsca;
out.grainradius=grainradius;
out.dust=dust;

fprintf('finished writing cube %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));
end

