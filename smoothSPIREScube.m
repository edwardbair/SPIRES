function out=smoothSPIREScube(nameprefix,vars,divisor,dtype,outloc,matdates,...
    grainradius_nPersist,mask,topofile,el_cutoff,fsca_thresh,cc,fice,...
    endconditions)
%function to smooth cube after running through SPIRES
% nameprefix - name prefix for outputs, e.g. Sierra
% vars - cell, variable list
% divisor - divisors for variables
% dtype - datatype for each variable
% outloc - output location, string
% matdates - datenum vector for image days
% grainradius_nPersist: min # of consecutive days needed with normal 
% grain sizes to be kept as snow, e.g. 7
% mask- logical mask w/ ones for areas to exclude 
% topofile- h5 file name from consolidateTopography, part of TopoHorizons
% el_cutoff, min elevation for snow, m - scalar, e.g. 1500
% fsca_thresh: min fsca cutoff, scalar e.g. 0.15
% cc - static canopy cover, single or doube, same size as watermask,
% 0-1 for viewable gap fraction correction
% fice - fraction of ice/neve, single or double, 0-1, mxn
% endcondition - string, end condition for splines for dust and grain size, 
% e.g. 'estimate' or 'periodic', see slmset.m

%output: struct out w/ fields
%fsca, grainradius, dust, and hdr (geographic info)

%1.34 hr/yr for h08v05, 2017
fprintf('reading %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

for i=1:length(matdates)
    dv=datevec(matdates(i));
    fname=fullfile(outloc,[nameprefix datestr(dv,'yyyymm') '.mat']);
    m=matfile(fname);
    if i==1
        for j=1:length(vars)
            out.(vars{j})=zeros([size(m.(vars{j}),1) size(m.(vars{j}),2) ...
                length(matdates)],'single');
        end
%         fsca=zeros([size(m.fsca_raw,1) size(m.fsca_raw,2) length(matdates)],'single');
%         fshade=zeros([size(m.fsca_raw,1) size(m.fsca_raw,2) length(matdates)],'single');
%         weights=zeros([size(m.fsca_raw,1) size(m.fsca_raw,2) length(matdates)],'single');
%         %note zeros due to uint storage
%         grainradius=zeros([size(m.fsca_raw,1) size(m.fsca_raw,2) length(matdates)],'single');
%         dust=zeros([size(m.fsca_raw,1) size(m.fsca_raw,2) length(matdates)],'single');    
    end
    ind=datenum(dv)-datenum([dv(1:2) 1])+1;%1st day of month

    %convert to single and scale
    for j=1:length(vars)
        tt=m.(vars{j})(:,:,ind)==intmax(dtype{j});
        v=single(m.(vars{j})(:,:,ind));
        v(tt)=NaN;
        v=v/divisor(j);
        out.(vars{j})(:,:,i)=v;
    end
    
    %fsca
    
%     fsca_t=single(m.fsca_raw(:,:,ind));
%     t=fsca_t==intmax('uint8');
%     fsca_t(t)=NaN;
%     fsca_t=fsca_t/100;
%     %fshade
%     fshade_t=single(m.fshade(:,:,ind));
%     t=fshade_t==intmax('uint8');
%     fshade_t(t)=NaN;
%     fshade_t=fshade_t/100;
%     %weights    
%     weights_t=single(m.weights(:,:,ind));
%     t=weights_t==intmax('uint8');
%     weights_t(t)=NaN;
%     weights_t=weights_t/100;
%     %grain radius
%     grainradius_t=single(m.grainradius(:,:,ind));
%     t=grainradius_t==intmax('uint16');
%     grainradius_t(t)=NaN;
%     %dust
%     dust_t=single(m.dust(:,:,ind));
%     t=dust_t==intmax('uint16');
%     dust_t(t)=NaN;
%     dust_t=dust_t/10;
%     
%     fsca(:,:,i)=fsca_t;
%     fshade(:,:,i)=fshade_t;
%     weights(:,:,i)=weights_t;
%     grainradius(:,:,i)=grainradius_t;
%     dust(:,:,i)=dust_t;
end

fprintf('finished reading %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

%store raw values before any adjustments
out.fsca_raw=out.fsca;

%fshade adj prior to smoothing
t=out.fshade<1;
out.fsca(t)=out.fsca(t)./(1-out.fshade(t));

%and viewable gap correction prior to smoothing

%enlarge cc by pixelsize

earthRadius = 6.371007181e+03;
orbitHeight = 705;
[ppl,ppt,~] = pixelSize(earthRadius,orbitHeight,1,out.sensorZ);

cc=cc.*ppl.*ppt;
cc(cc>1)=1;
cc(isnan(cc))=0;

out.fsca=out.fsca./(1-cc);
out.fsca(out.fsca>1)=1;
out.fsca(isnan(out.fsca))=0;

%elevation filter
[Z,hdr]=GetTopography(topofile,'elevation');
Zmask=Z < el_cutoff;
Zmask=repmat(Zmask,[1 1 length(matdates)]);
%masked area filter
bigmask=repmat(mask,[1 1 size(out.fsca,3)]);

out.fsca(Zmask | bigmask) = 0;

%create mask for cube where radius is > 50 & radius < 1190 for 7 or more days

gmask=snowPersistenceFilter(out.grainradius > 50 & out.grainradius < 1190,...
    grainradius_nPersist,1);

% set to NAN days that aren't in that mask but are not zero fsca
out.fsca(~gmask & ~(out.fsca==0))=NaN;

newweights=out.weights;
newweights(isnan(out.fsca))=0;

%fill in and smooth NaNs

fprintf('smoothing fsca %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

out.fsca=smoothDataCube(out.fsca,newweights,'mask',~mask,...
   'method','smoothingspline','SmoothingParam',0.02);

%get some small fsca values from smoothing - set to zero
out.fsca(out.fsca<fsca_thresh)=0;
out.fsca(bigmask)=NaN;

t0=out.fsca==0; %track zeros to prevent 0/0 = NaN

%need to rep matrix for logical operation below
fice=repmat(fice,[1 1 size(out.fsca,3)]);
% fice correction
fice(isnan(fice))=0;
out.fsca=out.fsca./(1-fice);

%set back to zero
out.fsca(t0)=0;
%fix high values
out.fsca(out.fsca>1)=1;

%fix values below thresh to ice values
t=out.fsca<fice;
out.fsca(t)=fice(t);

fprintf('finished smoothing fsca %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

fprintf('smoothing grain radius %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

%create mask of any fsca for interpolation
anyfsca=any(out.fsca,3);

newweights=out.weights;
newweights(isnan(out.fsca) | out.fsca==0)=0; %no snow, no weight
newweights(~gmask)=0;

dF=cat(3,zeros(size(out.fsca,1,2)),diff(out.fsca,1,3));
%send logical cube for decreasing fsca
fcube=dF<=0;
out.grainradius=smoothDataCube(out.grainradius,newweights,'mask',anyfsca,...
    'method','slm','monotonic','increasing','fcube',fcube,'knots',-3,...
    'endconditions',endconditions);

out.grainradius(out.fsca==0 | isnan(out.fsca))=NaN;

fprintf('finished smoothing grain radius %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

fprintf('smoothing dust %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

%compute grain sizes differences
dG=cat(3,zeros(size(out.grainradius,1,2)),diff(out.grainradius,1,3));
%send logical cube for increasing grain sizes
fcube=dG>=0;

out.dust=smoothDataCube(out.dust,newweights,'mask',anyfsca,...
    'method','slm','monotonic','increasing','fcube',fcube,'knots',-3,...
    'endconditions',endconditions);
out.dust(out.fsca==0 | isnan(out.fsca))=NaN;

fprintf('finished smoothing dust %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

fprintf('writing cube %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

% out.fsca_raw=fsca_raw; %store elevation filtered & masked 
% %(from input mask) fsca
out.matdates=matdates;
out.hdr=hdr;
% out.fsca=fsca;
% out.grainradius=grainradius;
% out.dust=dust;

fprintf('finished writing cube %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));
end

