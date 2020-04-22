function out=run_spires(R0,R,solarZ,Ffile,watermask,fsca_thresh,...
    pshade,dustmask,tolval,cc,hdr,red_b,swir_b)
% run LUT version of scagd for 4-D matrix R
% produces cube of: fsca, grain size (um), and dust concentration (by mass)
% input:
% R0 - background image (MxNxb). Recommend a sinle day with no clouds, 
%no snow and and low sensor zenith angle
% R - 4D cube of time-space smoothed reflectances (MxNxbxd) with
% dimensions: x,y,band,day
% solarZ: solar zenith angles (MxNxd) for R
% Ffile, location of griddedInterpolant object that produces reflectances
%for each band
% with inputs: grain radius (um), dust (ppmw), solar zenith angle (deg),
% band # (not necessarily in order of wavelength)
% watermask: logical mask (MxN), true for water
% fsca_thresh: min fsca cutoff, scalar e.g. 0.15
% pshade: shade spectra (bx1) or photometric scalar, e.g. 0
% dustmask: make where dust values can be retrieved, 0-1
% tolval: unique row tolerance value, i.e. 0.05 - bigger number goes faster
% as more pixels are grouped together
% cc - canopy cover (MxN), 0-1. No NaNs
% hdr - header struct w/ map info, see GetCoordinateInfo.m
% red_b - red band, e.g. 3 for MODIS and L8
% swir_b - SWIR band, e.g. 6 for MODIS and L8

%output:
%   out : struct w fields
%   fsca: MxNxd
%   grainradius: MxNxd
%   dust: MxNxd

sz=size(R);

if length(sz) == 3
   sz(4)=1; % singleton 4th dimenson
end

[X,Y]=pixcenters(hdr.RefMatrix,size(watermask),'makegrid');

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
cc=reshape(cc,[sz(1)*sz(2) 1]);
dm=reshape(dustmask,[sz(1)*sz(2) 1]);

for i=1:sz(4) %for each day
    thisR=squeeze(R(:,:,i));
    thissolarZ=squeeze(solarZ(:,i));
    
    NDSI=(thisR(:,red_b)-thisR(:,swir_b))./...
        (thisR(:,red_b)+thisR(:,swir_b));
    t=NDSI > 0  & ~watermask & ~isnan(thissolarZ);
    M=[round(thisR,2) round(R0,2) round(thissolarZ) round(cc,2) dm];
    
    %keep track of indices
    Rind=1:sz(3);
    R0ind=(Rind(end)+1):(Rind(end)+sz(3));
    sZind=R0ind(end)+1;
    ccind=sZind(end)+1;
    dmind=ccind(end)+1;
    
    M=M(t,:); % only values w/ > 0 NDSI and no water
    XM=X(t); % X coordinates for M
    YM=Y(t); % Y coordinates for M
    [c,im,~]=uniquetol(M,tolval,'ByRows',true,...
        'DataScale',1,'OutputAllIndices',true);
    t1=tic;
    temp=zeros(size(c,1),4); %fsca,shade,grain radius,dust
    
    
    parfor j=1:size(c,1) %solve for unique (w/ tol) rows
        pxR=c(j,Rind);
        pxR0=c(j,R0ind);
        sZ=c(j,sZind);
        thiscc=c(j,ccind);
        thisdm=c(j,dmind);
        o=speedyinvert(pxR,pxR0,sZ,Ffile,pshade,thisdm,[],thiscc);
        sol=o.x; %fsca,shade,grain radius,dust
        if sol(1) < 0.75 %set dust to NaN for low fsca_raw
           sol(4)=NaN;  
        end
        sol(1)=sol(1)/(1-sol(2));%normalize by fshade
        temp(j,:)=sol;
    end
    %make a copy of temp for use below
    temp2=temp;
    %re-solve for places w/ NaN dust 
    tt=~isnan(temp(:,4));
    if nnz(tt) >= 4 % if there are at least 4 solved dust values
        sdust=temp(tt,4); %solved dust values
        parfor j=1:size(c,1) 
            if isnan(temp(j,4)) % if there's no dust value, re-solve using
                %interpolated dust 
                pxR=c(j,Rind);
                pxR0=c(j,R0ind);
                sZ=c(j,sZind);
                thiscc=c(j,ccind);
                thisdm=c(j,dmind);
                idx=im{j}(1); %index to row of M correspnding to c
                dist=pdist2([XM(tt),YM(tt)],[XM(idx),YM(idx)]);
                %inverse distance weighted average
                w=1./dist;
                D=(sum(w.*sdust))./(sum(w));
                if D >= 0.1 %if snow is dirty
                    o=speedyinvert(pxR,pxR0,sZ,Ffile,pshade,thisdm,D,...
                        thiscc);
                    sol=o.x; %fsca,shade,grain radius,dust
                    sol(1)=sol(1)/(1-sol(2));%normalize by fshade
                    temp2(j,:)=sol;
                else
                    temp2(j,:)=[temp(j,1:3) 0];
                end  
            end
        end
    end
    %create a list of indices
    %and fill matrices corresponding to NDSI > 0 & ~water
    %can't use parfor for this
    repxx=zeros(size(M,1),4);
    for j=1:size(temp2,1) % the unique indices
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
    t2=toc(t1);
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
