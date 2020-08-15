function out=smoothSPIREScube(nameprefix,outloc,matdates,...
    windowSize,windowThresh,mingrainradius,maxgrainradius,mindust,maxdust,...
    mask,topofile,el_cutoff,fsca_thresh,cc,fice,b_R)
%function to smooth cube after running through SPIRES
% nameprefix - name prefix for outputs, e.g. Sierra
% outloc - output location, string
% matdates - datenum vector for image days
% windowSize - search window size for moving persistence filter, e.g. 45
% windowThresh - threshold number of days w/ fsca in windows to avoid being
% zeroed, e.g. 13
% mingrainradius - min believable grain radius, um, e.g. 75 um
% maxgrainradius - max believable grain radius, e.g. 1100 um
% mindust - min dust content, e.g. 0 um
% maxdust - max believable dust: max believable dust, e.g. 950 ppm
% mask- logical mask w/ ones for areas to exclude
% topofile - h5 file name from consolidateTopography, part of TopoHorizons
% el_cutoff - min elevation for snow, m - scalar, e.g. 1000
% fsca_thresh - min fsca cutoff, scalar e.g. 0.10
% cc - static canopy cover, single or doube, same size as mask,
% 0-1 for viewable gap fraction correction
% fice - fraction of ice/neve, single or double, 0-1, mxn
% b_R - b/R ratio for canopy cover, see GOvgf.m, e.g. 2.7

%output: struct out w/ fields
%fsca, grainradius, dust, and hdr (geographic info)

%0.75-1.5 hrfor Sierra Nevada with 60 cores
time1=tic;

fprintf('reading %s...%s\n',datestr(matdates(1)),datestr(matdates(end)));
%int vars
vars={'fsca','fshade','grainradius','dust','weights','sensorZ'};
divisor=[100 100 1 10 100 1];
dtype={'uint8','uint8','uint16','uint16','uint8','uint8'};

dv=datevec(matdates);
dv=dv(dv(:,3)==1,:);

for i=1:size(dv,1)
    fname=fullfile(outloc,[nameprefix datestr(dv(i,:),'yyyymm') '.mat']);
    m=matfile(fname);
    if i==1
        for j=1:length(vars)
            out.(vars{j})=zeros([size(m.(vars{j}),1) size(m.(vars{j}),2) ...
                length(matdates)],'single');
        end
    end
    doy_start=datenum(dv(i,:))-datenum(dv(1,:))+1;
    doy_end=doy_start+size(m.fsca,3)-1;

    %convert to single and scale
    for j=1:length(vars)
        tt=m.(vars{j})==intmax(dtype{j});
        v=single(m.(vars{j}));
        v(tt)=NaN;
        v=v/divisor(j);
        out.(vars{j})(:,:,doy_start:doy_end)=v;
    end
end

fprintf('finished reading %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

%store raw values before any adjustments
out.fsca_raw=out.fsca;

tmask=out.fsca>fsca_thresh & out.grainradius > mingrainradius;
tmask=movingPersist(tmask,windowSize,windowThresh);
out.fsca(~tmask)=0;

cc(isnan(cc))=0;
t=out.fsca==0;

%use GO model
cc_adj=1-GOvgf(cc,0,0,out.sensorZ,0,b_R);

fice(isnan(fice))=0;
fice=repmat(fice,[1 1 size(out.fsca,3)]);

%combine cc and fshade adjustment
out.fsca=out.fsca./(1-cc_adj-out.fshade-fice);
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

badg=out.grainradius<mingrainradius | out.grainradius>maxgrainradius | ...
out.dust > maxdust ;

%grain sizes too small or large to be trusted
out.grainradius(badg)=NaN;

[~,idx]=sort(out.grainradius,3,'descend','MissingPlacement','last');

%find indices of N ranked max values
N=3;
idx=idx(:,:,1:N);
%find latest occuring peak
idx=max(idx,[],3)+1;

%set everything after peak to that value
for i=1:size(idx,1)
    for j=1:size(idx,2)
        out.grainradius(i,j,idx(i,j):end)=out.grainradius(i,j,idx(i,j)-1);
    end 
end
    
%now fix dust values also before smoothing grain sizes
%set dust to zero for small grains
out.dust(badg)=NaN;

%do the same for dust
[~,idx]=sort(out.dust,3,'descend','MissingPlacement','last');
idx=idx(:,:,1:N);
idx=max(idx,[],3)+1;
for i=1:size(idx,1)
    for j=1:size(idx,2)
        out.dust(i,j,idx(i,j):end)=out.dust(i,j,idx(i,j)-1);
    end 
end

% use weights
newweights=out.weights;
newweights(isnan(out.grainradius) | out.fsca==0)=0;


out.grainradius=smoothDataCube(out.grainradius,newweights,'mask',anyfsca,...
   'method','smoothingspline','SmoothingParam',0.8);

fprintf('finished smoothing grain radius %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

out.grainradius(out.grainradius<mingrainradius)=mingrainradius;
out.grainradius(out.grainradius>maxgrainradius)=maxgrainradius;
out.grainradius(out.fsca==0)=NaN;

fprintf('smoothing dust %s...%s\n',datestr(matdates(1)),...
    datestr(matdates(end)));

%set dust to zero when smoothed grain radius is below some value
dust_rg_thresh=250; %um

out.dust(out.grainradius <= dust_rg_thresh) = mindust;
dv=datevec(matdates);
t=matdates<datenum([dv(end,1) 3 1]);
out.dust(:,:,t) = mindust;

out.dust=smoothDataCube(out.dust,newweights,'mask',anyfsca,...
     'method','smoothingspline','SmoothingParam',0.1);

%  out.dust=smoothDataCube(out.dust,newweights,'mask',anyfsca,...
%       'method','slm','knots',-15,'envelope','supremum');


% out.dust=smoothDataCube(out.dust,newweights,'mask',anyfsca,...
%     'method','slm','monotonic','increasing','fcube',out.dust>0,'knots',-2);
  
%clean up out of bounds splines
out.dust(out.dust>maxdust)=maxdust;
out.dust(out.dust<mindust)=mindust;
t=matdates<datenum([dv(end,1) 3 1]);
out.dust(:,:,t) = mindust;

out.dust(out.fsca==0)=NaN;

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

