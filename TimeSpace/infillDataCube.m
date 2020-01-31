function [ output_cube ] = infillDataCube( cube, varargin )
% [ output_cube ] = infillDataCube( cube [, mask] )
%
% Fills in NaN values in a 3D data cube
% Input
%   cube - raw cube, with NaN specifying missing values
% Optional input
%   mask (logical) - 2D portion of cube to be analyzed, if missing then whole cube
%
% Output
%   output_cube with NaN values filled in, original values not modified
%
% Jeff Dozier 2016-01-15
% modified by NB 2020-1-7
%% process required argument
p = inputParser;
defaultMask = false(0,0);
addRequired(p,'cube',@isnumeric)
addOptional(p,'mask',defaultMask,@islogical)
parse(p,cube,varargin{:})
cube = p.Results.cube;
% whole image if mask not specified
if isempty(p.Results.mask)
    mask = true(size(cube,1),size(cube,2));
    whole = true;
else
    mask = p.Results.mask;
    whole = false;
end

%% truncate cube to bounding box of mask
if ~whole
%     stats = regionprops(mask,'BoundingBox');
%     cc = floor(struct2array(stats.BoundingBox));
%     bb = [cc(2)-1 cc(1)-1 cc(2)+cc(4)+1 cc(1)+cc(3)+1];
%     bb(3:4) = min(bb(3:4),[size(cube,1) size(cube,2)]);
%     bb(1:2) = max(bb(1:2),[1 1]);
    [r,c]=find(mask);
    bb(1)=min(r);
    bb(3)=max(r);
    bb(2)=min(c);
    bb(4)=max(c);
    iCube = cube(bb(1):bb(3),bb(2):bb(4),:);
    iMask = mask(bb(1):bb(3),bb(2):bb(4));
else
    iCube = cube;
    iMask = mask;
end
% for passing to blockproc, concatenate the mask to the end of the cube
if whole
    limits = [nanmin(iCube(:)) nanmax(iCube(:))];
else
    rMask = repmat(iMask,[1 1 size(iCube,3)]);
    limits = [nanmin(iCube(rMask)) nanmax(iCube(rMask))];
end
iCube = cat(3,iCube,cast(iMask,'like',iCube));

% call blockproc
fun = @(block_struct) blockInterp(block_struct.data);
siz = bestblk([size(iCube,1) size(iCube,2)],200);
if isequal(siz,[size(iCube,1) size(iCube,2)])
    newCube = blockInterp(iCube); % don't use blockproc if whole cube at once
else
%     if ~isempty(invoke_parpool())
q=gcp('nocreate');
    if ~isempty(q)
        useParallel = true;
    else
        warning('not able to start parallel pool, running sequentially')
        useParallel = false;
    end
    newCube = blockproc(iCube,siz,fun,'BorderSize',[8 8],'useParallel',useParallel);
end

if whole
    output_cube = truncateLimits(newCube,limits);
else
    % replace bounding box part of original cube
    output_cube = cube;
    output_cube(bb(1):bb(3),bb(2):bb(4),:) = truncateLimits(newCube,limits);
end

end

function Z = blockInterp(XM)
% interpolate X, constrained by mask M

% separate the cube from the mask and replicate the mask to be the size of
% the cube
N = size(XM);
X = double(XM(:,:,1:N(3)-1));
M = repmat(logical(XM(:,:,N(3))),1,1,size(X,3));

% if all NaNs or no NaNs, just return
if ~any(isnan(X(M))) || all(isnan(X(M)))
    Z = cast(X,'like',XM);
    return
end
% if no variability, fill in and return
xu = unique(X);
if (length(xu)==2) % i.e., a value and NaN
    Z = cast(ones(size(X))*nanmax(xu),'like',XM);
    return
end

% points to interpolate from and to
[x,y,z] = meshgrid(1:N(2),1:N(1),1:N(3));
t = M & ~isnan(X);
ptsI = [x(t) y(t) z(t)];
V = X(t);
tn = M & isnan(X);
ptsQ = [x(tn) y(tn) z(tn)];
% add some columns to make a viable alpha shape
pt = min(ptsQ,[],1)-1;
pb = max(ptsQ,[],1)+1;
pt = max(pt,[1 1 1]);
pb = min(pb,[size(X,2) size(X,1) size(X,3)]);

% alpha shape for points to interpolate to, multiplying the default alpha
% radius by 1.5
pt = unique(cat(1,pt,ptsQ,pb),'rows');
shp = alphaShape(pt);
shp = alphaShape(pt,1.5*shp.Alpha);

% interpolating points in the alpha shape, and scattered Interp of those
try
    ti = inShape(shp,ptsI);
    S = scatteredInterpolant(ptsI(ti,:),V(ti));
    if nnz(ti)<10
        S = scatteredInterpolant(ptsI,V);
    end
catch
    S = scatteredInterpolant(ptsI,V);
end

% interpolate and put in output
Z = X;
zhat = S(ptsQ);
try
    Z(tn) = zhat;
catch
    S = scatteredInterpolant(ptsI,V);
    zhat = S(ptsQ);
    Z(tn) = zhat;
end

Z = cast(Z,'like',XM);

end