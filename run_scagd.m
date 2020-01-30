function out=run_scagd_modis(R0,R,solarZ,Ffile,watermask,fsca_thresh,...
    pshade,dust_thresh,tolval)
% run LUT version of scagd for 4-D matrix R
% produces cube of: fsca, grain size (um), and dust concentration (by mass)
% input:
% R0 - background image (MxNxb). Recommend using time-spaced smoothed
% cube from a month with minimum fsca and clouds, like August or September,
% then taking minimum of reflectance for each band (b)
% R - 4D cube of time-space smoothed reflectances (MxNxbxd) with
% dimensions: x,y,band,day
% solarZ: solar zenith angles (MxNxd) for R
% Ffile, location of griddedInterpolant object that produces reflectances
%for each band
% with inputs: grain radius, dust, solar zenith angle, band
% watermask: logical mask, true for water
% fsca_thresh: min fsca cutoff, scalar e.g. 0.15
% pshade: shade spectra (bx1); reflectances
% dust_tresh: threshold cutoff to return dust values, e.g. 0.85
% tolval: unique row tolerance value, i.e. 0.05 - bigger number goes faster
% as more pixels are grouped together

%output:
%   out : struct w fields
%   fsca: MxNxd
%   grainradius: MxNxd
%   dust: MxNxd

sz=size(R);

if length(sz) == 3
   sz(4)=1; % singleton 4th dimenson
end

[X,Y]=meshgrid(1:sz(1),1:sz(2));


fsca=zeros([sz(1)*sz(2) sz(4)]);
grainradius=NaN([sz(1)*sz(2) sz(4)]);
dust=NaN([sz(1)*sz(2) sz(4)]);

solarZ=reshape(double(solarZ),[sz(1)*sz(2) sz(4)]);
R=reshape(double(R),[sz(1)*sz(2) sz(3) sz(4)]);
R0=reshape(double(R0),[sz(1)*sz(2) sz(3)]);
watermask=reshape(watermask,[sz(1)*sz(2) 1]);
shade=zeros(size(fsca));
X=reshape(X,[sz(1)*sz(2) 1]);
Y=reshape(Y,[sz(1)*sz(2) 1]);

red_b=3;
swir_b=6;

for i=1:sz(4) %for each day
    thisR=squeeze(R(:,:,i));
    thissolarZ=squeeze(solarZ(:,i));
    
    NDSI=(thisR(:,red_b)-thisR(:,swir_b))./...
        (thisR(:,red_b)+thisR(:,swir_b));
    t=NDSI > 0  & ~watermask & ~isnan(thissolarZ);
    M=[round(thisR,2) round(R0,2) round(thissolarZ)];
    M=M(t,:); % only values w/ > 0 NDSI and no water
    XM=X(t); % X coordinates for M
    YM=Y(t); % Y coordinates for M
    [c,im,~]=uniquetol(M,tolval,'ByRows',true,...
        'DataScale',1,'OutputAllIndices',true);
    tic;
    temp=zeros(length(c),4); %fsca,shade,grain radius,dust
    parfor j=1:length(c) %solve for unique (w/ tol) rows
        pxR=c(j,1:7);
        pxR0=c(j,8:14);
        sZ=c(j,15);
        o=speedyinvert(pxR,pxR0,sZ,Ffile,pshade,dust_thresh,[]);
        sol=o.x; %fsca,shade,grain radius,dust
        sol(1)=sol(1)/(1-sol(2));%normalize by fshade
        temp(j,:)=sol;
    end
    %make a copy of temp for use below
    temp2=temp;
    %median_dustval=median(temp(:,4),'omitnan');
    %re-solve for places w/ NaN dust 
    tt=~isnan(temp(:,4));
    if nnz(tt) >= 4 % if there are at least 4 solved dust values
        sdust=temp(tt,4); %solved dust values
        parfor j=1:length(c) 
            if isnan(temp(j,4)) % if there's no dust value, re-solve using
                %interpolated dust 
                pxR=c(j,1:7);
                pxR0=c(j,8:14);
                sZ=c(j,15);
                idx=im{j}(1); %index to row of M correspnding to c
                [~,idx_s]=pdist2([XM(tt),YM(tt)],[XM(idx),YM(idx)],...
                    'euclidean','Smallest',4);
                D=mean(sdust(idx_s)); % mean of closest 4 dust values
                o=speedyinvert(pxR,pxR0,sZ,Ffile,pshade,dust_thresh,D);
                sol=o.x; %fsca,shade,grain radius,dust
                sol(1)=sol(1)/(1-sol(2));%normalize by fshade
                temp2(j,:)=sol;
            end
        end
    end
    %create a list of indices
    %and fill matrices corresponding to NDSI > 0 & ~water
    %can't use parfor for this
    repxx=zeros(length(M),4);
    for j=1:length(temp2) % the unique indices
        idx=im{j}; %indices for each unique val
        repxx(idx,1)=temp2(j,1); %fsca
        %don't forget shade is elem 2
        repxx(idx,2)=temp2(j,2); %shade
        repxx(idx,3)=temp2(j,3); %grain radius
        repxx(idx,4)=temp2(j,4); %dust
    end
    %now fill out all pixels
    fsca(t,i)=repxx(:,1);
    shade(t,i)=repxx(:,2);
    grainradius(t,i)=repxx(:,3);
    dust(t,i)=repxx(:,4);
    t2=toc;
    fprintf('done w/ day %i in %g min\n',i,t2/60);
end

fsca=reshape(fsca,[sz(1) sz(2) sz(4)]);
grainradius=reshape(grainradius,[sz(1) sz(2) sz(4)]);
dust=reshape(dust,[sz(1) sz(2) sz(4)]);
shade=reshape(shade,[sz(1) sz(2) sz(4)]);

fsca(fsca<fsca_thresh)=0;
grainradius(fsca==0)=NaN;
dust(fsca==0)=NaN;
shade(fsca==0)=NaN;

out.fsca=fsca;
out.grainradius=grainradius;
out.dust=dust;
out.shade=shade;
