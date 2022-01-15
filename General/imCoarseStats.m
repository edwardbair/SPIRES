function [Cbar,Cstddev,Ccount]=imCoarseStats(Cmap,Csize,Fmap,F)
% [Cbar,Cstddev,Ccount]=imCoarseStats(Cmap,Csize,Fmap,F)
% calculates mean and std dev for image with coarse spatial resolution,
% based on image with fine spatial resolution
% [Cbar,Cstddev,Ccount]=imCoarseStats(Cmap,Csize,Fmap,F)
% input
%   Cmap - referencing matrix for coarse-resolution image
%   Csize - size (row,col) for coarse-resolution image
%   Fmap - referencing matrix for fine-resolution image
%   F - the fine-resolution image
% output
%   Cbar - coarse-resolution image with mean values from F
%   Cstddev - coarse-resolution image with std dev values from F
%   Ccount - image with count of fine pixels for each cell
%
%   Note: values of Cbar and Cstddev are set to NaN if they are outside
%   the bounding box of the F image.
%
% output matrices same as input, but must be floating point

assert(isfloat(F),'input F matrix must be single or double')
Cbar=NaN(Csize,'like',F);
Cstddev=NaN(Csize,'like',F);
Ccount=zeros(Csize,'like',F);
%
Crow=Csize(1);
Ccol=Csize(2);
[Frow,Fcol]=size(F);
[Cx,Cy]=pixcenters(Cmap,Csize,'makegrid');
[Fx,Fy]=pixcenters(Fmap,size(F)); % don't use makegrid here because of memory
%
% threshold size, squared
thresh=max([abs(Cmap(2,1)) abs(Cmap(1,2))])/2;
thresh=(sqrt(2)*thresh)^2;
%
% pixel coordinates of F for each coarse-resolution pixel
[Fprow,Fpcol]=map2pix(Fmap,Cx,Cy);
Fprow=int32(round(Fprow));
Fpcol=int32(round(Fpcol));
%
% window size in the fine-resolution image
multiplier=round(max([abs(Cmap(2,1)) abs(Cmap(1,2))])/...
    min([abs(Fmap(2,1)) abs(Fmap(1,2))]));
% loop thru the coarse-resolution image
for j=1:Crow
    for k=1:Ccol
        % window in F that is close to pixel j,k in coarse image
        minrow=max([1 Fprow(j,k)-multiplier]);
        maxrow=min([Fprow(j,k)+multiplier Frow]);
        mincol=max([1 Fpcol(j,k)-multiplier]);
        maxcol=min([Fpcol(j,k)+multiplier Fcol]);
        %
        % subset of F
        WF=F(minrow:maxrow,mincol:maxcol);
        [Wy,Wx]=ndgrid(Fy(minrow:maxrow),Fx(mincol:maxcol));
        %
        % which values in WF are within threshold?
        t=(Cx(j,k)-Wx).^2+(Cy(j,k)-Wy).^2 < thresh;
        Cbar(j,k)=mean(WF(t),'omitnan');
        Cstddev(j,k)=std(WF(t),'omitnan');
        Ccount(j,k)=nnz(~isnan(WF(t)));
    end
end