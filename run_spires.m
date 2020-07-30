function out=run_spires(R0,R,solarZ,Ffile,mask,shade,tolval,...
    red_b,swir_b,solarZthresh)
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
% mask: logical mask (MxN), true for areas to NOT process
% shade endmember, scalar or vector, length of # bands
% tolval: unique row tolerance value, i.e. 0.05 - bigger number goes faster
% as more pixels are grouped together
% red_b - red band, e.g. 3 for MODIS and L8
% swir_b - SWIR band, e.g. 6 for MODIS and L8
% solarZthresh - max solar zenith angle for dust estimates; 
% otherwise assumed clean, e.g. 40 deg

%output:
%   out : struct w fields
%   fsca: MxNxd
%   fshade: MxNxd
%   grainradius: MxNxd
%   dust: MxNxd

outvars={'fsca','fshade','grainradius','dust'};

sz=size(R);

if length(sz) == 3
   sz(4)=1; % singleton 4th dimenson
end

for i=1:length(outvars)
    out.(outvars{i})=NaN([sz(1)*sz(2) sz(4)]);
end

solarZ=reshape(double(solarZ),[sz(1)*sz(2) sz(4)]);
R=reshape(double(R),[sz(1)*sz(2) sz(3) sz(4)]);
R0=reshape(double(R0),[sz(1)*sz(2) sz(3)]);
mask=reshape(mask,[sz(1)*sz(2) 1]);

for i=1:sz(4) %for each day
    thisR=squeeze(R(:,:,i));
    thissolarZ=squeeze(solarZ(:,i));
    
    NDSI=(thisR(:,red_b)-thisR(:,swir_b))./...
         (thisR(:,red_b)+thisR(:,swir_b));
    t=NDSI > -0.5  & ~mask & ~isnan(thissolarZ) & all(~isnan(thisR),2);
  
    M=[round(thisR,2) round(R0,2) round(thissolarZ)];    
    %keep track of indices
    Rind=1:sz(3);
    R0ind=(Rind(end)+1):(Rind(end)+sz(3));
    sZind=R0ind(end)+1;
    
    M=M(t,:);
    [c,im,~]=uniquetol(M,tolval,'ByRows',true,'DataScale',1,...
        'OutputAllIndices',true);
    im1=zeros(size(im));
    for k=1:length(im)
       im1(k)=im{k}(1); 
    end
        
    t1=tic;

    temp=NaN(size(c,1),length(outvars));
    
    parfor j=1:size(c,1) %solve for unique (w/ tol) rows

        pxR=c(j,Rind);
        pxR0=c(j,R0ind);
        sZ=c(j,sZind);
        if sZ > solarZthresh %assume clean
            o=speedyinvert(pxR,pxR0,sZ,Ffile,shade,0,1,[],[]);
        else %assume dirty
            o=speedyinvert(pxR,pxR0,sZ,Ffile,shade,1,1,[],[]);
        end
        temp(j,:)=o.x;
    end

    repxx=NaN(size(M,1),length(outvars));
    for j=1:size(temp,1) % the unique indices
        idx=im{j}; %indices for each unique val
        for k=1:length(outvars)
            repxx(idx,k)=temp(j,k); %fsca,fshade,grain size,dust
        end
    end
    %now fill out all pixels
    for j=1:length(outvars)
        out.(outvars{j})(t,i)=repxx(:,j);
        tt=NDSI <= -0.5;
        if strcmp(outvars{j}(1),'f') %set fsca/fshade to zero
            out.(outvars{j})(tt,i)=0;
        else
            out.(outvars{j})(tt,i)=NaN;%grain & dust to NaN
        end
    end
    t2=toc(t1);
    fprintf('done w/ day %i in %g min\n',i,t2/60);
end

for i=1:length(outvars)
    out.(outvars{i})=reshape(out.(outvars{i}),[sz(1) sz(2) sz(4)]);
end

end