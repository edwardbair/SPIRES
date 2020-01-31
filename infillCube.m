function [ newcube] = infillCube( cube,weight,varargin)
% [ newcube] = infillCube( cube, weight,varargin )
%fill in the NaNs in a time-space data cube, but not where weight==0 or
%snow_fraction==0
%
% input
%   cube - 3D space-time data cube
%   weight - confidence weights
% optional input
%   snow_fraction mask (i.e. snow_fraction>0)
%
% output
%   newcube - input cube with NaNs filled
%
% Jeff Dozier, 2015-08-16
%
% to parallelize the infilling across multiple cores, run with block sizes
% that include a buffer
linDim = 500;
bufferSize = 50;
matrixThreshold = linDim^2;
matrixSize = [size(cube,1) size(cube,2)];
if matrixSize(1)*matrixSize(2)>matrixThreshold
    blksiz = bestblk(matrixSize,linDim-1);
else
    blksiz = matrixSize;
end

limits = [nanmin(cube(:)) nanmax(cube(:))];
if nargin>2 %optional arg is a mask
    cube = setMask(cube,varargin{1},-Inf); % set to infinity to indicate we don't interpolate
end

% inpaint the NaNs and create next version of spacetime cube
% to ameliorate edge effects, run with a buffer
assert(isequal(size(cube),size(weight)),['sizes of cubes ('...
    num2str(size(cube)) ') and weights (' num2str(size(weight)) ') must be same'])
if isequal(blksiz,matrixSize)
    Z = fill_NaNs(cube,weight,blksiz);
else
    Z = fill_NaNs(cube,weight,blksiz,bufferSize);
end
newcube = cast(truncateLimits(Z,limits),'like',cube);
end

function ns = fill_NaNs(X,W,siz,varargin)
% inpaint the NaNs, using blockproc

% can pass just one array to block_struct, so combine weights with X
W = cast(W,'like',X);
XW = cat(3,X,W);
fun = @(block_struct) blockScatInterp(block_struct.data);

if nargin>3
    bufferSize = varargin{1};
    ns = blockproc(XW,siz,fun,'UseParallel',true,'BorderSize',...
        [bufferSize bufferSize]);
else
    ns = blockproc(XW,siz,fun,'UseParallel',true);
end
ns = cast(ns,'like',X);
end

function Z = blockScatInterp(XW)
% uses MATLAB scatteredInterpolant

% separate weights from values
N = size(XW,3);
assert(mod(N,2)==0,'something wrong, size(XW,3) should be even')
N = N/2;
X = XW(:,:,1:N);
W = XW(:,:,N+1:end);

% if all NaNs or no NaNs, just return
if all(isnan(X(:)) | isinf(X(:))) ||...
        ~any(isnan(X(:)) | isinf(X(:))) ||...
        all(W(:)==0)
    X(isinf(X)) = NaN;
    Z = X;
    return
end
% if all zeros, just return
if nansum(X(~isinf(X(:))))==0
    Z = zeros(size(X),'like',X);
    Z(isinf(X) | W==0) = NaN;
    return;
end

% coordinates
[x,y,z] = meshgrid(1:size(X,2),1:size(X,1),1:size(X,3));
N = size(X);
V = double(X(:));
W = double(W(:));
x = x(:);
y = y(:);
z = z(:);
toFill = isnan(V) & W>0;
toLeave = W==0 | isinf(V);
toBound = V>0 & ~toLeave;
stayZero = V==0 & ~toLeave;

% identify the alpha shape for the points to consider and fill
if any(toBound(:))
    shp = alphaShape(x(toBound),y(toBound),z(toBound));
    if ~isempty(shp) && size(shp.Points,1)>=4 % threshold, 4 points
        try
            tPoss = inShape(shp,x,y,z);
            toFill = tPoss & toFill;
            % use the zeros inside the alpha shape
            if any(toFill(:))
                toUse = tPoss & ~(toFill | toLeave);
                F = scatteredInterpolant(x(toUse),y(toUse),z(toUse),V(toUse));
                V(toFill) = F(x(toFill),y(toFill),z(toFill));
            end
        catch
            toLeave = W==0 | isinf(V); % dummy, sort of
        end
    elseif any(toFill)
        toUse = ~(toFill | toLeave);
        if any(toUse)
            V(toFill) = nanmedian(V(toUse));
        else
            V(toFill) = 0;
        end
    end
end
if any(toLeave)
    V(toLeave) = NaN;
end
if any(stayZero)
    V(stayZero) = 0;
end
Z = cast(reshape(V,N),'like',X);
end