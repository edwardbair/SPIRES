function out=smoothSPIREScube(nameprefix,outloc,matdates,...
    nPersistDry,nPersistSnow,mingrainradius,maxgrainradius,mindust,maxdust,...
    mask,topofile,el_cutoff,fsca_thresh,cc,fice,...
    endconditions)
%function to smooth cube after running through SPIRES
% nameprefix - name prefix for outputs, e.g. Sierra
% outloc - output location, string
% matdates - datenum vector for image days
% nPersistDry: min # of consecutive days to trust a dry land retrieval,
% e.g. 4
% nPersistSnow: min # of consectuve days to trust a snow covered
% retreival,e.g. 8, remember that these persistence filters are performed 
% on cloud (and other gap) filled estimates
% mingrainradius: min believable grain radius, um, e.g. 100 um
% maxgrainradius: max believable grain radius, e.g. 1100 um
% mindust: min dust content, e.g. 12 um
% maxdust: max believable dust: max believable dust, e.g. 950 ppm
% mask- logical mask w/ ones for areas to exclude
% topofile- h5 file name from consolidateTopography, part of TopoHorizons
% el_cutoff, min elevation for snow, m - scalar, e.g. 1000
% fsca_thresh: min fsca cutoff, scalar e.g. 0.10
% cc - static canopy cover, single or doube, same size as watermask,
% 0-1 for viewable gap fraction correction
% fice - fraction of ice/neve, single or double, 0-1, mxn
% endcondition - string, end condition for splines for dust and grain size, 
% e.g. 'estimate' or 'periodic', see slmset.m

%output: struct out w/ fields
%fsca, grainradius, dust, and hdr (geographic info)

%1.34 hr/yr for h08v05, 2017
time1=tic;

fprintf('reading %s...%s\n',datestr(matdates(1)),datestr(matdates(end)));
%int vars
vars={'fsca','fshade','grainradius','dust','weights','sensorZ'};
divisor=[100 100 1 10 100 1];
dtype={'uint8','uint8','uint16','uint16','uint8','uint8'};

for i=1:length(matdates)
    dv=datevec(matdates(i));
    fname=fullfile(outloc,[nameprefix datestr(dv,'yyyymm') '.mat']);
    m=matfile(fname);
    if i==1
        for j=1:length(vars)
            out.(vars{j})=zeros([size(m.(vars{j}),1) size(m.(vars{j}),2) ...
                length(matdates)],'single');
        end
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
end

fprintf('finished reading %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

%store raw values before any adjustments
out.fsca_raw=out.fsca;

%create trust mask for zero fsca values, skipped or missing is NaN
tmask=out.fsca<fsca_thresh & out.weights>0;
dmask=snowPersistenceFilter(tmask,nPersistDry,1);

out.fsca(dmask)=0;

%create trust mask for non-zero snow values
tmask=out.fsca>fsca_thresh;
smask=snowPersistenceFilter(tmask,nPersistSnow,1);

out.fsca(~smask & ~dmask)=NaN;

%fshade adj 
% t=out.fshade<1;
% out.fsca(t)=out.fsca(t)./(1-out.fshade(t));

%and viewable gap correction 

%enlarge cc by pixelsize
% earthRadius = 6.371007181e+03;
% orbitHeight = 705;
% [ppl,ppt,~] = pixelSize(earthRadius,orbitHeight,1,out.sensorZ);
% 
% cc=cc.*ppl.*ppt;
% cc(cc>1)=1;

cc(isnan(cc))=0;
t=out.fsca==0;
% t1=out.fsca>fsca_thresh & cc>0;
% out.fsca(t1)=1;
% out.fsca=out.fsca./(1-cc);
% out.fsca(out.fsca>1)=1; %includes cc==1 case of Inf, which is any pos #/0
% out.fsca(t)=0;

%White et al 2005, Figure 1ab, https://doi.org/10.1080/01431160500080626
a=-0.67;
b=0.089;
cc_adj1=a.*cc+b;

a=-0.79;
b=0.074;
cc_adj2=a.*cc+b;

%negative indicates underestimate
cc_adj=cc-min(cat(3,cc_adj1,cc_adj2),[],3);

%combine cc and fshade adjustment
out.fsca=out.fsca./(1-cc_adj-out.fshade);
out.fsca(out.fsca>1 | out.fsca<0)=1;
out.fsca(t)=0;

%elevation filter
[Z,hdr]=GetTopography(topofile,'elevation');
Zmask=Z < el_cutoff;
Zmask=repmat(Zmask,[1 1 length(matdates)]);

%masked area filter
bigmask=repmat(mask,[1 1 size(out.fsca,3)]);

out.fsca(Zmask | bigmask) = 0;

newweights=out.weights;
newweights(isnan(out.fsca))=0;

%fill in and smooth NaNs

fprintf('smoothing fsca %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

out.fsca=smoothDataCube(out.fsca,newweights,'mask',~mask,...
   'method','smoothingspline','SmoothingParam',0.1);

%get some small fsca values from smoothing - set to zero
out.fsca(out.fsca<fsca_thresh)=0;
out.fsca(bigmask)=NaN;

% t0=out.fsca==0; %track zeros to prevent 0/0 = NaN

%need to rep matrix for logical operation below
fice(isnan(fice))=0;
fice=repmat(fice,[1 1 size(out.fsca,3)]);
% fice correction

% out.fsca=out.fsca./(1-fice);

%set back to zero
% out.fsca(t0)=0;
%fix high values
% out.fsca(out.fsca>1)=1;

%fix values below thresh to ice values
t=out.fsca<fice;
out.fsca(t)=fice(t);
out.fsca(out.fsca<fsca_thresh)=0;

fprintf('finished smoothing fsca %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

fprintf('smoothing grain radius %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

%create mask of any fsca for interpolation
anyfsca=any(out.fsca,3);

%send logical cube for decreasing fsca
dF=cat(3,zeros(size(out.fsca,1,2)),diff(out.fsca,1,3));
fcube=dF<0;

%grain sizes too small or large to be trusted
%bad grain sizes
badg=out.grainradius<mingrainradius | out.grainradius>maxgrainradius;
out.grainradius(badg)=NaN;

% use weights
newweights=out.weights;
newweights(isnan(out.grainradius) | out.fsca==0)=0;

out.grainradius=smoothDataCube(out.grainradius,newweights,'mask',anyfsca,...
    'method','slm','monotonic','increasing','fcube',fcube,'knots',-2,...
    'endconditions',endconditions);

out.grainradius(out.grainradius<mingrainradius)=mingrainradius;
out.grainradius(out.grainradius>maxgrainradius)=maxgrainradius;

out.grainradius(out.fsca==0 | isnan(out.fsca))=NaN;

fprintf('finished smoothing grain radius %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

fprintf('smoothing dust %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

% assume bad grain sizes are bad dust values
badd=out.dust>maxdust;
out.dust(badg | badd)=NaN;

%median filter for dust
% use weights
newweights=out.weights;
newweights(isnan(out.dust) | out.fsca==0)=0;

% t=~isnan(out.dust);
% out.dust(t)=movmedian(out.dust(t),movfiltlength,3);

out.dust=smoothDataCube(out.dust,newweights,'mask',anyfsca,...
    'method','slm','monotonic','increasing','fcube',fcube,'knots',-2,...
    'endconditions',endconditions);

out.dust(out.dust>maxdust)=maxdust;
out.dust(out.dust<mindust)=mindust;
out.dust(out.fsca==0 | isnan(out.fsca))=NaN;

fprintf('finished smoothing dust %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

%write out h5 cubes
fname=fullfile(outloc,[nameprefix datestr(matdates(end),'yyyy') '.h5']);
if exist(fname,'file')
    delete(fname); 
end

out.matdates=matdates;
out.hdr=hdr;

fprintf('writing cubes %s...%s\n',datestr(matdates(1)),...
datestr(matdates(end)));

%output variables
outvars={'fsca_raw','fsca','grainradius','dust'};
outnames={'raw_snow_fraction','snow_fraction','grain_size','dust'};
outdtype={'uint8','uint8','uint16','uint16'};
outdivisors=[100 100 1 10];

for i=1:length(outvars)   
    member=outnames{i};
    Value=out.(outvars{i});
    dS.(member).divisor=outdivisors(i);
    dS.(member).dataType=outdtype{i};
    dS.(member).maxVal=max(Value(:));
    dS.(member).FillValue=intmax(dS.(member).dataType);
    writeh5stcubes(fname,dS,out.hdr,out.matdates,member,Value);
end

time2=toc(time1);
fprintf('completed in %5.2f hr\n',time2/60/60);

end

