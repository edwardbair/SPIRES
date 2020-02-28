function [ output_cube ] = smoothDataCube(cube, weight,varargin)
%
% Time-smooths data cube either with smoothing spline or discrete cosine transform
% Input
%   cube - raw cube
%   weight - corresponding weights
% Optional input - name-value pairs
%   'mask' (logical) - 2D portion of cube to be analyzed, if missing then whole cube
%   'method' - either 'smoothingspline', 'slm' (shape language toolbox), 'smoothn' (default, Damien Garcia's
%       smoothn from the MATLAB file exchange)
%   Note: smoothn is a simultaneous 3D smoothing (i.e. spatial smearing),
%   while smoothingspline only smooths in the time dimension
%   'bordersize' (m x n), default 8x8, see "doc blockproc.m"
%   for slm method only:
%   'monotonic', string either 'increasing' or 'decreasing, see slmset.m
%   'fcube' (m x n), difference datacube to follow for monotonic 
%    increasing/decreasing, see slmset.m
%   'knots',number of knots for spline, see slmset.m
%   'endconditions', e.g. 'estimate' (default), see slmset.m
%
% Output
%   output_cube - smoothed values, NaN outside mask if specified
%
% Jeff Dozier 2016-01-23
% modified by NB 2020-02-17

%% process required argument
p = inputParser;
defaultMask = false(0,0);
defaultMethod = 'smoothn';
defaultBorderSize = [8 8];
defaultBlockSize = bestblk([size(cube,1) size(cube,2)],200);
defaultSmoothingParam= [];
defaultfcube = [];
defaultMonotonic = 'increasing';
defaultKnots = 6;
defaultEndConditions = 'estimate';

addRequired(p,'cube',@isnumeric)
addRequired(p,'weight',@isnumeric)
addParameter(p,'mask',defaultMask,@islogical)
addParameter(p,'method',defaultMethod,@ischar)
addParameter(p,'bordersize',defaultBorderSize,@isnumeric);
addParameter(p,'blocksize',defaultBlockSize,@isnumeric);
addParameter(p,'smoothingparam',defaultSmoothingParam,@isnumeric);
addParameter(p,'monotonic',defaultMonotonic,@ischar);
addParameter(p,'fcube',defaultfcube,@(x)ismatrix(x)||islogical(x));
addParameter(p,'knots',defaultKnots,@isnumeric);
addParameter(p,'endconditions',defaultEndConditions,@ischar);

parse(p,cube,weight,varargin{:})
cube = p.Results.cube;
weight = p.Results.weight;
sp=p.Results.smoothingparam;
monotonic=p.Results.monotonic;
fcube=p.Results.fcube;
knots=p.Results.knots;
endconditions=p.Results.endconditions;

assert(isequal(size(cube),size(weight)),'size of cube and weights must be equal')
% whole image if mask not specified
if isempty(p.Results.mask)
    mask = true(size(cube,1),size(cube,2));
    %whole = true;
else
    mask = p.Results.mask;
    %whole = false;
end

method = p.Results.method;
bordersize = p.Results.bordersize;
blocksize = p.Results.blocksize;

N = size(cube);

option.TolZ = 1.e-3;

switch method
    case 'smoothn_sequential'        
        %reshape and transpose cube to put time sequence next to each other in memory
        % (MATLAB is column-major order)
        cube = reshape(cube,N(1)*N(2),N(3))';
        weight = reshape(weight,N(1)*N(2),N(3))';
        mask = reshape(mask,N(1)*N(2),1)';
        % boundaries of the mask
        fcol = find(mask,1,'first');
        if isempty(fcol) % if the mask is all false
            output_cube = reshape(cube',N(1),N(2),N(3));
            return
        end
        lcol = find(mask,1,'last');
        iCube = double(cube(:,fcol:lcol));
        weight = double(weight(:,fcol:lcol));
        mask = mask(fcol:lcol);
        limits = [nanmin(iCube(:)) nanmax(iCube(:))];
        
        % by column
        sCube = nan(size(iCube));
        weight(isnan(weight) | weight==0) = .01; % make sure we have enough weights
        k=gcp('nocreate');
        if ~isempty(k)
            parfor c=1:size(iCube,2)
                if mask(c)
                    y = iCube(:,c);
                        t=~isnan(y);
                    if nansum(y) == 0 || nnz(t)<2
                        sCube(:,c) = zeros(size(y));
                    else
                        % smooth
                        tooSmall = 'MATLAB:smoothn:SLowerBound';
                        tooLarge = 'MATLAB:smoothn:SUpperBound';
                        maxIter = 'MATLAB:smoothn:MaxIter';
                        warning('off',tooSmall);
                        warning('off',tooLarge);
                        warning('off',maxIter);
                        [Z,~,flag] = smoothn(y,weight(:,c),'robust',option);
                    if ~flag
                        x=1:length(y);
                        x=x';
                        ex=isnan(y);
                        F = fit(x,y,'smoothingspline','weights',weight(:,c),...
                            'exclude',ex);
                        %Z = smooth(y,'rlowess');
                        Z = F(x);
                    end
                        Z = cast(Z,'like',y);
                        sCube(:,c) = Z;
                    end
                end
            end
        else
            error('start a parpool up');
        end
        
        % back into original
        cube(:,fcol:lcol) = cast(truncateLimits(sCube,limits),'like',cube);
        
        output_cube = reshape(cube',N(1),N(2),N(3));
        
    case 'smoothingspline'
        %% reshape and transpose cube to put time sequence next to each other in memory
        % (MATLAB is column-major order)
        cube = reshape(cube,N(1)*N(2),N(3))';
        weight = reshape(weight,N(1)*N(2),N(3))';
        mask = reshape(mask,N(1)*N(2),1)';
        % boundaries of the mask
        fcol = find(mask,1,'first');
        if isempty(fcol) % if the mask is all false
            output_cube = reshape(cube',N(1),N(2),N(3));
            return
        end
        lcol = find(mask,1,'last');
        iCube = double(cube(:,fcol:lcol));
        weight = double(weight(:,fcol:lcol));
        mask = mask(fcol:lcol);
        limits = [nanmin(iCube(:)) nanmax(iCube(:))];
        
        % by column
        sCube = nan(size(iCube));
        weight(isnan(weight) | weight==0) = .01; % make sure we have enough weights
        x = (1:size(iCube,1))';
        %         if ~isempty(invoke_parpool())
        k=gcp('nocreate');
        if ~isempty(k)
            parfor c=1:size(iCube,2)
                if mask(c)
                    y = iCube(:,c);
                        t=~isnan(y);
                    if nansum(y) == 0 || nnz(t)<2
                        sCube(:,c) = zeros(size(y));
                    else
                        t=isnan(y);
                        F = fit(x,y,'smoothingspline','weights',weight(:,c),...
                            'exclude',t,'SmoothingParam',sp);
                        sCube(:,c) = F(x);
                    end
                end
            end
        else
            error('start a parpool');
        end
        
        % back into original
        cube(:,fcol:lcol) = cast(truncateLimits(sCube,limits),'like',cube);
        
        output_cube = reshape(cube',N(1),N(2),N(3));
        
    case 'slm'
        % reshape and transpose cube to put time sequence next to each other in memory
        % (MATLAB is column-major order)
        monotonicflag=true;
        if isempty(fcube)
            monotonicflag=false;
        end
        cube = reshape(cube,N(1)*N(2),N(3))';
        weight = reshape(weight,N(1)*N(2),N(3))';
        mask = reshape(mask,N(1)*N(2),1)';
        % boundaries of the mask
        fcol = find(mask,1,'first');
        if isempty(fcol) % if the mask is all false
            output_cube = reshape(cube',N(1),N(2),N(3));
            return
        end
        lcol = find(mask,1,'last');
        iCube = double(cube(:,fcol:lcol));
        if monotonicflag
            fcube = reshape(fcube,N(1)*N(2),N(3))';
            FCube = double(fcube(:,fcol:lcol));
        else
            FCube = zeros(size(iCube)); %needed for parfor
        end
        weight = double(weight(:,fcol:lcol));
        mask = mask(fcol:lcol);
        limits = [nanmin(iCube(:)) nanmax(iCube(:))];
        
        % by column
        sCube = nan(size(iCube));
        weight(isnan(weight) | weight==0) = .01; % make sure we have enough weights
        x = (1:size(iCube,1))';
        k=gcp('nocreate');
        if ~isempty(k)
            parfor c=1:size(iCube,2)
                if mask(c)
                    y = iCube(:,c);
                        t=~isnan(y);
                    if nansum(y) == 0 || nnz(t)<2
                        sCube(:,c) = zeros(size(y));
                    else
                        if monotonicflag
                            f=FCube(:,c);
                            %find contiguous regions
                            [start,finish]=contiguous(f);
                            %use those to specify where spline should be
                            %monotonically increasing/decreasing
                            F=slmengine(x,double(y),monotonic,...
                                [start finish],'weights',weight(:,c),...
                                'knots',knots,'endconditions',endconditions);
                            sCube(:,c) = slmeval(x,F);
                        else
                            F=slmengine(x,double(y),'weights',weight(:,c),...
                                'knots',knots,'endconditions',endconditions);
                            sCube(:,c) = slmeval(x,F);
                        end
                    end
                end
            end
        else
            error('start a parpool');
        end
        
        % back into original
        cube(:,fcol:lcol) = cast(truncateLimits(sCube,limits),'like',cube);
        
        output_cube = reshape(cube',N(1),N(2),N(3));
        
    case 'smoothn'
        iCube = cube;
        iMask = mask;
        % for passing to blockproc, concatenate the weights and the mask to the end of the cube
        limits = [nanmin(iCube(iMask)) nanmax(iCube(iMask))];
        iCube = cat(3,iCube,weight,cast(iMask,'like',iCube));
        
        % call blockproc
        fun = @(block_struct) blockSmooth(block_struct.data);
        siz = blocksize;
%         siz = bestblk([size(iCube,1) size(iCube,2)],200);
        
        if isequal(siz,[size(iCube,1) size(iCube,2)]) % don't use blockproc if just 1 block
            output_cube = blockSmooth(iCube);
        else
            k=gcp('nocreate');
            if ~isempty(k)
                %~isempty(invoke_parpool())
                useParallel = true;
            else
                warning('not able to start parallel pool, running sequentially')
                useParallel = false;
            end
            newCube = blockproc(iCube,siz,fun,'BorderSize',bordersize,...
                'useParallel',useParallel);
            
            % replace bounding box part of original cube
            output_cube = truncateLimits(newCube,limits);
        end
        warning('on','all')
        
    otherwise
        error('method %s not recognized',method)
end
end

function Z = blockSmooth(XM)
%% smooth X with weight W constrained by mask M
% turn off the warnings in smoothn, as the blockSmooth function
% uses the exit flag
tooSmall = 'MATLAB:smoothn:SLowerBound';
tooLarge = 'MATLAB:smoothn:SUpperBound';
maxIter = 'MATLAB:smoothn:MaxIter';
warning('off',tooSmall);
warning('off',tooLarge);
warning('off',maxIter);
% separate the cube from the weights and mask and replicate the mask to be the size of
% the cube
N = size(XM);
assert(mod(N(3),2)==1,'3rd dimension of cube size must be odd, N=%d,%d,%d',...
    N(1),N(2),N(3))
endData = (N(3)-1)/2;
X = double(XM(:,:,1:endData));
W = double(XM(:,:,endData+1:N(3)-1));
M = repmat(logical(XM(:,:,N(3))),1,1,size(X,3));

% if all NaNs, just return
if all(isnan(X(M)))
    Z = cast(X,'like',XM);
    return
end
% if no variability, fill in and return
xu = unique(X);
if (length(xu)==2) % i.e., a value and NaN
    Z = cast(ones(size(X))*nanmax(xu),'like',XM);
    return
end

% smooth
option.TolZ = 1.e-2;
% if bad exit flag, just do simple smoothing
[Z,~,flag] = smoothn(X,W,'robust',option);
if ~flag
    Z = smooth3(X,'gaussian');
    Z(isnan(Z)) = X(isnan(Z));
end
Z = cast(Z,'like',XM);

end