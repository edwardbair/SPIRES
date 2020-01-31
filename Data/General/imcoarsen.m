function [ Rmap, Cbar ] = imcoarsen( Cmap, Fmap, F, varargin )
% [ Rmap, Cbar ] = imcoarsen( Cmap, Fmap, F [, optional input] )
%imcoarsen coarsens fine-scale image to match coarse-resolution image
% (for images that have a referencing matrix)
%
% Input
%   Cmap - referencing matrix for coarse-resolution image
%   Fmap - referencing matrix for fine-resolution image
%   F - the fine resolution image
% Optional input, one or both, in either order
%   'extra' - flag for extra coarsening
%       if omitted, coarsens to resolution of coarse image
%       if present, coarsens 1 step beyond, then refines back
%   method - resampling method, either 'nearest', 'bilinear', or 'bicubic'
%       (default is 'bilinear')
%
% Output
%   Rmap - referencing matrix for the coarsened image (not the same as Cmap
%       because the coarsened image covers only the area of F, not the
%       original area associated with Cmap)
%   Cbar - the coarsened image (data type same as F)

% size of the coarse image that corresponds to the area of F, based on
% coordinates of F corners (cw from NW)
sf = size(F);
[x,y] = pix2map(Fmap,[1 1 sf(1) sf(1)],[1 sf(2) sf(2) 1]);
Fx = [min(x) max(x)];
Fy = [max(y) min(y)];
[r,c] = map2pix(Cmap,Fx,Fy);
Csize = [round(1+max(r)-min(r)) round(1+max(c)-min(c))];

% preserve range
fmin = min(reshape(F,numel(F),1));
fmax = max(reshape(F,numel(F),1));

% optional arguments
optargin=size(varargin,2);
assert( nargin-optargin >= 3,...
    'need 3 required arguments: Cmap, Fmap, F')
% defaults
rstop = 2;
method = 'bilinear';
if optargin>0
    assert(optargin<3,'just 1 or 2 optional arguments permitted')
    for k=1:optargin
        if strcmpi(varargin{k},'extra')
            rstop = 1; % extra coarsening
        else
            method = varargin{k}; % imresize will return error if not recognized
        end
    end
end

% iteratively resize F
% (if extra, resizes to larger pixel than Cmap, then contracts back)
Cbar = F;
while max(size(Cbar)./Csize)>rstop
    Cnew = imresize(Cbar,1/2,method);
    if any(isnan(Cnew(:))) && ~strcmpi(method,'nearest')
        Cnear = imresize(Cbar,1/2,'nearest');
        t = isnan(Cnew);
        Cnew(t) = Cnear(t);
    end
    Cbar = Cnew;
end
% final resize, and preserve range
Cnew = imresize(Cbar,Csize,method);
if any(isnan(Cnew(:))) && ~strcmpi(method,'nearest')
    Cnear = imresize(Cbar,Csize,'nearest');
    t = isnan(Cnew);
    Cnew(t) = Cnear(t);
end
Cbar = Cnew;

t = Cbar<fmin;
Cbar(t) = fmin;
t = Cbar>fmax;
Cbar(t) = fmax;

% make sure pixel size is correct
pixelsize = (size(F)./Csize).*[Fmap(2,1) Fmap(1,2)]; % check size OK
reldiff = (pixelsize-[Cmap(2,1) Cmap(1,2)])./[Cmap(2,1) Cmap(1,2)];
assert(max(abs(reldiff))<0.05,...
    'relative pixel differences %g %g too different - debug',reldiff)
% output referencing matrix
x11 = Fx(1)+Cmap(2,1)/2;
y11 = Fy(1)+Cmap(1,2)/2;
Rmap = makerefmat(x11,y11,Cmap(2,1),Cmap(1,2));

end