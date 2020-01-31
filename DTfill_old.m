function [ Z,varargout ] = DTfill( X,mask,varargin )
% [ Z ] = DTfill( X,mask [,'nDays',value,'threshold',value] )
% fill NaNs in masked area of 3D cube using Delaunay triangulation
%
% input
%   X - input cube
%   mask - area mask of x,y dimensions of X, or 3D mask of size(X)
% optional input, name-value pairs
%   'nDays' - only consider pixels that have consecutive nDays of X value
%       (default empty)
%   'threshold' - value below which set to zero (default 0.0)
%
% output
%   Z - filled cube with NaN outside area of mask filled
% optional output
%   W - adjusted weights to keep track of what's been filled

defaultThreshold = 0;
defaultDays = [];
p = inputParser;
addRequired(p,'X',@isnumeric)
addRequired(p,'mask',@islogical)
addParameter(p,'threshold',defaultThreshold,@isnumeric)
addParameter(p,'nDays',defaultDays,@isnumeric)
parse(p,X,mask,varargin{:})
threshold = p.Results.threshold;
nDays = p.Results.nDays;

N = size(X);
[x,y,z] = meshgrid(1:N(2),1:N(1),1:N(3));
% replicate the 2D mask to 3D if not already
if ~isequal(size(mask),N)
    assert(isequal(size(mask),N(1:2)),'mask must be size of x-y dimensions of X')
    mask = repmat(mask,1,1,N(3));
end
tnan = ~mask | isnan(X);
t = X>0 & mask;
X(~t) = 0;
X(tnan) = NaN;
if ~isempty(nDays)
    X = snowPersistenceFilter(X,nDays,threshold);
end

if nnz(X>0)==0
    Z = zeros(size(X),'like',X);
    Z(~mask) = NaN;
    if nargout>1
        W = ones(size(X));
        W(~mask) = 0;
        varargout{1} = W;
    end
    return
end

% interpolate from the non-zero values inside the mask, and their
% coordinates, but restricted to pixels with nDays duration if specified
t = X>0 & mask;
ptsI = [x(t) y(t) z(t)];
V = X(t);
% initialize output, keeping all NaNs
Z = X;
% Delaunay triangulation of the points to interpolate from
DT = delaunayTriangulation(ptsI);
% alpha shape of the points
shp = alphaShape(ptsI);
% points to interpolate to - the NaN values inside the alpha shape
tn = mask & isnan(X);
ptsQ = [x(tn) y(tn) z(tn)];
tf = inShape(shp,ptsQ);
ptsQ = ptsQ(tf,:);
% nearest neighbors to unknown points (distance to DT is always smaller)
IDT = nearestNeighbor(DT,ptsQ);
% values at nearest neighbors
for k=1:length(IDT)
    Z(ptsQ(k,2),ptsQ(k,1),ptsQ(k,3)) = V(IDT(k));
end
% any remaining NaN values inside the mask are zero
Z(isnan(Z)) = 0;
Z(~mask) = NaN;
% run again thru the persistence filter
if ~isempty(nDays)
    Z = snowPersistenceFilter(Z,nDays,threshold);
end

% reset weights to account for points filled, inverse to distance to
% nearest neighbor
if nargout>1
    wt = zeros(size(ptsQ,1),1);
    for k=1:length(wt)
        wt(k) = 1./sqrt(sum((ptsQ(k,:)-DT.Points(IDT(k),:)).^2));
    end
    W = ones(size(Z));
    W(~mask) = 0;
    for k=1:size(ptsQ,1)
        W(ptsQ(k,2),ptsQ(k,1),ptsQ(k,3)) = wt(k);
    end
    % final check in case we divided by zero (should not have)
    W(isinf(W)) = 1;
    varargout{1} = W;
end

end