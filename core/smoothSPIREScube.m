function out=smoothSPIREScube(nameprefix,outloc,matdates,...
    windowSize,windowThresh,mingrainradius,maxgrainradius,mindust,maxdust,...
    mask,topofile,el_cutoff,fsca_thresh,cc,fice,b_R,dust_rg_thresh,...
    fixpeak,Nd)
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
% dust_rg_thresh, min grain radius for dust, e.g. 400 um
% maxflag
% fixpeak - boolean, true fixes grain and dust values at after peak at 
% peak value. Avoids physically impossible retrievals such as shrinking fsca 
% and grain size due to increasing SWIR reflectance
% Nd - number of days from length(matdates) to stop fixing peak, ignored if
% fixpeak is false
%output: struct out w/ fields
%fsca, grainradius, dust, and hdr (geographic info)

%0.75-1.5 hrfor Sierra Nevada with 60 cores
time1=tic;

h5name=fullfile(outloc,[nameprefix datestr(matdates(end),'yyyy') '.h5']);

if exist(h5name,'file')
    fprintf('%s exists, skipping\n',h5name);
else
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
    
    %run binary fsca mask through temporal filter
    tmask=out.fsca>fsca_thresh & out.grainradius > mingrainradius;
    tmask=movingPersist(tmask,windowSize,windowThresh);
    
    %create 2 smoothed versions: fsca (adjusted for cc,ice,shade,
    % elevation cutoff,watermask, fsca_min)
    %and fsca_raw (no cc,ice adj, or shade adj), but elevation cutoff, watermask, &
    %fsca_min applied)
    out.fsca(~tmask)=0;
    out.fsca_raw(~tmask)=0;
    
    cc(isnan(cc))=0;
    t=out.fsca==0;
    
    %use GO model
    cc_adj=1-GOvgf(cc,0,0,out.sensorZ,0,b_R);
    
    fice(isnan(fice))=0;
    fice=repmat(fice,[1 1 size(out.fsca,3)]);
    
    %combine cc and fshade adjustment
    out.fsca=out.fsca./(1-cc_adj-out.fshade-fice);
    out.fsca(out.fsca>1 | out.fsca<0)=1;
    %fix 0/0
    out.fsca(t)=0;
    
    %elevation filter
    [Z,hdr]=GetTopography(topofile,'elevation');
    Zmask=Z < el_cutoff;
    Zmask=repmat(Zmask,[1 1 length(matdates)]);
    
    %masked area filter
    bigmask=repmat(mask,[1 1 size(out.fsca,3)]);
    
    out.fsca(Zmask | bigmask) = 0;
    out.fsca_raw(Zmask | bigmask) = 0;
   
    newweights=out.weights; 
    newweights(isnan(out.fsca))=0;
    
    %fill in and smooth NaNs
    
    fprintf('smoothing fsca,fsca_raw,fshade %s...%s\n',datestr(matdates(1)),...
        datestr(matdates(end)));
    
    %smooth fully adj fsca
    out.fsca=smoothDataCube(out.fsca,newweights,'mask',~mask,...
        'method','smoothingspline','SmoothingParam',0.1);
    
    %smooth fsca_raw
    out.fsca_raw=smoothDataCube(out.fsca_raw,newweights,'mask',~mask,...
        'method','smoothingspline','SmoothingParam',0.1);

    %smooth fshade
    out.fshade=smoothDataCube(out.fshade,newweights,'mask',~mask,...
        'method','smoothingspline','SmoothingParam',0.1);

    %get some small fsca values from smoothing - set to zero
    out.fsca(out.fsca<fsca_thresh)=0;
    out.fsca(bigmask)=NaN;
    
    %same for fsca_raw
    out.fsca_raw(out.fsca_raw<fsca_thresh)=0;
    out.fsca_raw(bigmask)=NaN;
    
    %same for fshade
    out.fshade(out.fsca_raw<fsca_thresh)=0;
    out.fshade(bigmask)=NaN;

    %fix values below thresh to ice values
    t=out.fsca<fice;
    out.fsca(t)=fice(t);
    out.fsca(out.fsca<fsca_thresh)=0;
    
    fprintf('finished smoothing fsca,fsca_raw,fshade %s...%s\n',datestr(matdates(1)),...
        datestr(matdates(end)));
    
    fprintf('smoothing grain radius %s...%s\n',datestr(matdates(1)),...
        datestr(matdates(end)));
    
    %create mask of any fsca for interpolation
    anyfsca=any(out.fsca,3);
    
    badg=out.grainradius<mingrainradius | out.grainradius>maxgrainradius | ...
        out.dust > maxdust ;
    
    %grain sizes too small or large to be trusted
    out.grainradius(badg)=NaN;
    out.dust(badg)=NaN;
    
    %grain sizes after melt out
    out.dust(out.fsca==0)=NaN;
    
    %find max grain size
    [~,idx]=sort(out.grainradius,3,'descend','MissingPlacement','last');
    
    %find indices of N ranked max values
    N=3;
    idx=idx(:,:,1:N);
    %find latest occuring peak
    idx=max(idx,[],3);
    %Nd=15; %N days to taper
    if fixpeak % set values after peak grain radius to peak   
    for i=1:size(idx,1)
        for j=1:size(idx,2)
            %get rid of spikes & drops
            rgvec=squeeze(out.grainradius(i,j,:));
            rgvec=hampel(rgvec,2,2);

            %day w/ max grain size
            maxDay=idx(i,j);
            %final day in the cube
            %lastDay=size(out.grainradius,3);
            
            %default endday is the final day in the cube
            endDay=size(out.grainradius,3)-Nd;
 
            %set those days to (near) max grain size
            ind=maxDay:endDay;
            
            maxrg=rgvec(maxDay);
            rgvec(ind)=maxrg;
            
            %set final values to temp values
            out.grainradius(i,j,:)=rgvec;

            %dust: set dust for all days prior to 0 if below grain thresh
            dustvec=squeeze(out.dust(i,j,:));
            
            %all days with small grain sizes
            tt=rgvec<=dust_rg_thresh; 
            
            %all days prior to max grain size
            ttstart=false(size(tt));
            ttstart(1:maxDay)=true;
            
            %all days prior to max grain size with small grains
            tt=ttstart & tt;
            
            %set dust on those days to zero
            dustvec(tt)=0;
            
            %use dust value from max rg day
            dval=out.dust(i,j,maxDay);
            
            %set dust after those days to value on maxday
            dustvec(ind)=dval;
            
            out.dust(i,j,:)=dustvec;
        end
    end
    
    else %don't fix values after peak grain size, 
        %assume zero dust for small grains
        out.dust(out.grainradius<dust_rg_thresh)=0; 
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

    if fixpeak
        %dust is always zero on day 1
        % force periodicity at 1 and end (usually 365/366)
        out.dust(:,:,1)=zeros(size(anyfsca));
        %slow smoothing from day 2 to Nd-end 
        out.dust(:,:,2:(end-Nd))=smoothDataCube(out.dust(:,:,2:(end-Nd)),...
            newweights(:,:,2:(end-Nd)),'mask',anyfsca,...
            'method','smoothingspline','SmoothingParam',0.1);
        %then taper to zero from Nd-end to end
        out.dust(:,:,(end-Nd):end)=taperDataCube(out.dust(:,:,(end-Nd):end),...
            zeros(size(anyfsca)),anyfsca);
    else
        out.dust=smoothDataCube(out.dust,newweights,'mask',anyfsca,...
            'method','smoothingspline','SmoothingParam',0.1);
    end
    %clean up out of bounds splines
    out.dust(out.dust>maxdust)=maxdust;
    out.dust(out.dust<mindust)=mindust;  
    out.dust(out.fsca==0)=NaN;
    
    fprintf('finished smoothing dust %s...%s\n',datestr(matdates(1)),...
        datestr(matdates(end)));
    
    %write out h5 cubes
    
    out.matdates=matdates;
    out.hdr=hdr;
    
    fprintf('writing cubes %s...%s\n',datestr(matdates(1)),...
        datestr(matdates(end)));
    
    %output variables
    outvars={'fsca_raw','fsca','fshade','grainradius','dust'};
    outnames={'raw_snow_fraction','snow_fraction','shade_fraction','grain_size','dust'};
    outdtype={'uint8','uint8','uint8','uint16','uint16'};
    outdivisors=[100 100 100 1 10];
    
    for i=1:length(outvars)
        member=outnames{i};
        Value=out.(outvars{i});
        dS.(member).divisor=outdivisors(i);
        dS.(member).dataType=outdtype{i};
        dS.(member).maxVal=max(Value(:));
        dS.(member).FillValue=intmax(dS.(member).dataType);
        writeh5stcubes(h5name,dS,out.hdr,out.matdates,member,Value);
    end
    
    time2=toc(time1);
    fprintf('completed in %5.2f hr\n',time2/60/60);
end
end