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
lockname=fullfile(outloc,[nameprefix datestr(matdates(end),'yyyy') '.h5lock']);

if exist(h5name,'file')==2
    fprintf('%s exists, skipping\n',h5name);
elseif exist(lockname,'file')==2
    fprintf('%s locked, skipping\n',lockname);
else
    fid=fopen(lockname,'w');
    fclose(fid);
    %delete lockname on cleanup
    cleanup=onCleanup(@()CleanupFun(lockname));

    fprintf('reading %s...%s\n',datestr(matdates(1)),datestr(matdates(end)));
    %int vars
    vars={'fsca','fshade','grainradius','dust','weights','sensorZ'};
    divisor=[100 100 1 10 100 1];
    dtype={'uint8','uint8','uint16','uint16','uint8','uint8'};

    dv=datevec(matdates);
    dv=dv(dv(:,3)==1,:);
    %check that full set of matdates exists
    for i=1:size(dv,1)
        fname=fullfile(outloc,[nameprefix datestr(dv(i,:),'yyyymm') '.mat']);
        %if fname doesn't exist,  delete lock, throw error
        if exist(fname,'file')==0
            delete(lockname);
            error('matfile %s doesnt exist\n',fname);
        end
    end
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

    fprintf('smoothing grain radius and dust %s...%s\n',datestr(matdates(1)),...
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
    %don't set out.grainradius to nan where fsca==0 until later
    %this helps maintain high grain size values

    % create new weights for grain size and dust
    newweights=out.weights;
    newweights(isnan(out.grainradius) | out.fsca==0)=0;

    %nearest neighbor interp. for areas with fsca, but all nan grain size or dust 
    anygrain=any(out.grainradius,3);
    anydust=any(out.dust,3);
    
    %fsca but no grain size values
    t=~anygrain & anyfsca;
    [r,c]=find(t);
    for i=1:length(r)
        X=[c(i),r(i)];
        [r2,c2]=find(anygrain);
        Y=[c2,r2];
        d=pdist2(X,Y);
        [~,idm]=min(d);
        rgfill=squeeze(out.grainradius(r2(idm),c2(idm),:));
        out.grainradius(r(i),c(i),:)=rgfill;
    end
    
    %fsca but no dust values
    t=~anydust & anyfsca;
    [r,c]=find(t);
    for i=1:length(r)
        X=[c(i),r(i)];
        [r2,c2]=find(anydust);
        Y=[c2,r2];
        d=pdist2(X,Y);
        [~,idm]=min(d);
        dustfill=squeeze(out.dust(r2(idm),c2(idm),:));
        out.dust(r(i),c(i),:)=dustfill;
    end
    


    if fixpeak % set values after peak grain radius to peak
        N=size(out.grainradius);
        %reshape to days x pixels
        grainradius=reshape(out.grainradius,N(1)*N(2),N(3))';
        dust=reshape(out.dust,N(1)*N(2),N(3))';
        fsca=reshape(out.fsca,N(1)*N(2),N(3))';
        weights=reshape(newweights,N(1)*N(2),N(3))';

        parfor i=1:size(grainradius,2)
            fscavec=squeeze(fsca(:,i));
            rgvec=squeeze(grainradius(:,i));
            weightsvec=squeeze(weights(:,i));

            %get rid of spikes & drops
            rgvec=hampel(rgvec,2,2);
            if max(fscavec)==0 %skip pixels w/ no snow
                continue;
            end
            t=fscavec>0;
            %last day with snow, or end of cube
            meltOutday=min(find(t,1,'last'),length(rgvec));
            %peak fixing cannot last more than N days from final peak
            %make a temp copy of rgvec
            rgvec_t=rgvec;
            maxFixedDays=40;%days
            %set all days prior to meltOutday-maxFixedDays to nan
            rgvec_t(1:(meltOutday-maxFixedDays))=nan;
            %find the max day & val
            [gs_maxVal,maxDay]=max(rgvec_t,[],'omitnan');
            %but if its nan, use the latest measured 
            if isnan(gs_maxVal)
                maxDay=find(~isnan(rgvec),1,'last');
            end
           
            endDay=length(rgvec)-Nd;

            %set those days to (near) max grain size
            ind=maxDay:endDay;

            maxrg=rgvec(maxDay);
            rgvec(ind)=maxrg;

            %smooth up to maxDay
            ids=1:maxDay-1;

            %set 1st day to min, may be set to nan later, but helps w/
            %keeping spline in check
            rgvec(1)=mingrainradius;
            weightsvec(1)=1;

            rgvec(ids)=smoothVector(ids',rgvec(ids),weightsvec(ids),0.8);
            %taper vector to min value
            grainradius(:,i)=taperVector(rgvec,Nd,mingrainradius);

            %dust: set dust for all days prior to 0 if below grain thresh
            dustvec=squeeze(dust(:,i));

            %all days with small grain sizes
            tt=rgvec<=dust_rg_thresh;

            %all days prior to max grain size
            ttstart=false(size(tt));
            ttstart(1:maxDay-1)=true;

            %all days prior to max grain size with small grains
            tt=ttstart & tt;

            %set dust on those days to zero
            dustvec(tt)=0;

            %use dust value from max rg day
            dval=dustvec(maxDay);
            
            %but if that value is nan, use max dust val
            if isnan(dval)
                dval=max(dustvec,[],'omitnan');
            end

            %set dust after those days to value on maxday
            dustvec(ind)=dval;
 
            %set dust to zero on day 1
            dustvec(1)=0;
            weightsvec(1)=1;
            %smooth up until maxday
            dustvec(ids)=smoothVector(ids',dustvec(ids),...
                weightsvec(ids),0.1);
            dust(:,i)=taperVector(dustvec,Nd,mindust);
        end
        %put back into cube
        out.fsca = reshape(fsca',N(1),N(2),N(3));
        out.grainradius = reshape(grainradius',N(1),N(2),N(3));
        out.dust = reshape(dust',N(1),N(2),N(3));
        
    else %don't fix values after peak grain size
        out.grainradius=smoothDataCube(out.grainradius,newweights,'mask',anyfsca,...
            'method','smoothingspline','SmoothingParam',0.8);
        %assume zero dust for small grains
        out.dust(out.grainradius<dust_rg_thresh)=0;
        out.dust=smoothDataCube(out.dust,newweights,'mask',anyfsca,...
            'method','smoothingspline','SmoothingParam',0.1);
    end

    fprintf('finished smoothing grain radius and dust %s...%s\n',datestr(matdates(1)),...
        datestr(matdates(end)));

    out.grainradius(out.grainradius<mingrainradius)=mingrainradius;
    out.grainradius(out.grainradius>maxgrainradius)=maxgrainradius;
    out.grainradius(out.fsca==0)=NaN;

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

    %create h5 cube in tmp then move to avoid network h5 write issues
    h5tmpname=fullfile(tempdir,[nameprefix datestr(matdates(end),'yyyy') '.h5']);
    for i=1:length(outvars)
        member=outnames{i};
        Value=out.(outvars{i});
        dS.(member).divisor=outdivisors(i);
        dS.(member).dataType=outdtype{i};
        dS.(member).maxVal=max(Value(:));
        dS.(member).FillValue=intmax(dS.(member).dataType);
        writeh5stcubes(h5tmpname,dS,out.hdr,out.matdates,member,Value);
    end
    system(['mv ' h5tmpname ' ' h5name]);
    delete(lockname);
    time2=toc(time1);
    fprintf('completed in %5.2f hr\n',time2/60/60);
end
end

function CleanupFun(lockname)
if exist(lockname,'file')==2
    fprintf('cleaning up %s\n',lockname)
    delete(lockname)
end
end