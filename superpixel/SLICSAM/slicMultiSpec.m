function [L,N,d] = slicMultiSpec(X, k, m, nItr,SAMconst)
% SLIC Simple Linear Iterative Clustering SuperPixels
%
% ADAPTED TO CLUSTER USING ALL BANDS OF MULTISPECTRAL SATELLITES AND
% TEXTURE FEATURES AS WELL
%
% Implementation of Achanta, Shaji, Smith, Lucchi, Fua and Susstrunk's
% SLIC Superpixels
%
% Usage:   [l, Am, Sp, d] = slic(im, k, m, seRadius, colopt, mw)
%
% Arguments:  X - Multi-Band Image to be segmented.
%              k - Number of desired superpixels. Note that this is nominal
%                  the actual number of superpixels generated will generally
%                  be a bit larger, especially if parameter m is small.
%              m - Weighting factor between colour and spatial
%                  differences. Values from about 5 to 40 are useful.  Use a
%                  large value to enforce superpixels with more regular and
%                  smoother shapes. Try a value of 10 to start with.
%       seRadius - Regions morphologically smaller than this are merged with
%                  adjacent regions. Try a value of 1 or 1.5.  Use 0 to
%                  disable.
%         colopt - String 'mean' or 'median' indicating how the cluster
%                  colour centre should be computed. Defaults to 'mean'
%             mw - Optional median filtering window size.  Image compression
%                  can result in noticeable artifacts in the a*b* components
%                  of the image.  Median filtering can reduce this. mw can be
%                  a single value in which case the same median filtering is
%                  applied to each L* a* and b* components.  Alternatively it
%                  can be a 2-vector where mw(1) specifies the median
%                  filtering window to be applied to L* and mw(2) is the
%                  median filtering window to be applied to a* and b*.
%
% Returns:     l - Labeled image of superpixels. Labels range from 1 to k.
%             Am - Adjacency matrix of segments.  Am(i, j) indicates whether
%                  segments labeled i and j are connected/adjacent
%             Sp - Superpixel attribute structure array with fields:
%                   L  - Mean L* value
%                   a  - Mean a* value
%                   b  - Mean b* value
%                   r  - Mean row value
%                   c  - Mean column value
%                   stdL  - Standard deviation of L*
%                   stda  - Standard deviation of a*
%                   stdb  - Standard deviation of b*
%                   N - Number of pixels
%                   edges - List of edge numbers that bound each
%                           superpixel. This field is allocated, but not set,
%                           by SLIC. Use SPEDGES for this.
%              d - Distance image giving the distance each pixel is from its
%                  associated superpixel centre.
%
% It is suggested that use of this function is followed by SPDBSCAN to perform a
% DBSCAN clustering of superpixels.  This results in a simple and fast
% segmentation of an image.
%
% Minor variations from the original algorithm as defined in Achanta et al's
% paper:
%
% - SuperPixel centres are initialised on a hexagonal grid rather than a square
%   one. This results in a segmentation that will be nominally 6-connected
%   which hopefully facilitates any subsequent post-processing that seeks to
%   merge superpixels.
% - Initial cluster positions are not shifted to point of lowest gradient
%   within a 3x3 neighbourhood because this will be rendered irrelevant the
%   first time cluster centres are updated.
%
% Reference: R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua and
% S. Susstrunk. "SLIC Superpixels Compared to State-of-the-Art Superpixel
% Methods"  PAMI. Vol 34 No 11.  November 2012. pp 2274-2281.
%
% See also: SPDBSCAN, MCLEANUPREGIONS, REGIONADJACENCY, DRAWREGIONBOUNDARIES, RGB2LAB

% Copyright (c) 2013 Peter Kovesi
% www.peterkovesi.com/matlabfns/
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% Feb  2013
% July 2013 Super pixel attributes returned as a structure array

% Note that most of the computation time is not in the clustering, but rather
% in the region cleanup process.

%UPDATES -tstillinger
% im - can be a multispectral image with numerous bands - not converted to
% L*a*b colorspace, kept as RTOA.
% colopt - add ability to cluster using spectral angle.
% can add a texture dimension to wieght and add to superpixel caluclation
%   imTex - texture image of gabor features
%   t - wieght to give texture features in classification
% IMAGE CANNOT HAVE FILL PIXLES RIGHT NOW
%removed edge wieght distance option
%working on making run with multispectral images with numerious fill pixels around
%boarder (OLI Level 1 WRS-2 aquistions)
%Changes teh distance function to work with SAM - sqrt was used wron gin
%prior implementation

%add ability to flag nan fill pixels and remove then add as a single labled
%region at end aslo - parfor anywhere???

%set all fill pixels as cluster 1, then skip for the rest of the
%computations

%fill pixels are L == -1

if ~exist('nItr','var')   || isempty(nItr),     nItr = 10;     end

[rows, cols, specBands] = size(X);

fillPxls = isnan(X);
fill2D = fillPxls(:,:,1);
X(fillPxls)=-9999;%ADDD SEE IF WORKS

% Nominal spacing between grid elements assuming hexagonal grid
S = sqrt(rows*cols / (k * sqrt(3)/2));

% Get nodes per row allowing a half column margin at one end that alternates
% from row to row
nodeCols = round(cols/S - 0.5);
% Given an integer number of nodes per row recompute S
S = cols/(nodeCols + 0.5);

% Get number of rows of nodes allowing 0.5 row margin top and bottom
nodeRows = round(rows/(sqrt(3)/2*S));
vSpacing = rows/nodeRows;

% Recompute k
k = nodeRows * nodeCols;

% Allocate memory and initialise clusters, labels and distances.

% Cluster centre data  1:bands is mean spectra and texture bands
% bands+1 is x (col), bands+2 is y (row) of centre, bands+3 is No of pixels

bands = specBands;
ccData = bands+2+1;  %-NEW METHOD DONT NEED LAST ROW if parfor used
C = zeros(ccData,k);

L = -ones(rows, cols);   % Pixel labels.
d = inf(rows, cols);     % Pixel distances from cluster centres.

% Initialise clusters on a hexagonal grid
kk = 1;
r = vSpacing/2;

for ri = 1:nodeRows %IMPROVE - SET UP SO FILL PIXELS ALONE ARE IN BORDER SUPERPIXELS
    % Following code alternates the starting column for each row of grid
    % points to obtain a hexagonal pattern. Note S and vSpacing are kept
    % as doubles to prevent errors accumulating across the grid.
    if mod(ri,2), c = S/2; else, c = S;  end
    
    for ci = 1:nodeCols
        cc = round(c); rr = round(r);
        
        if fill2D(rr,cc)
            c = c+S;
            
            continue
        else
            
            if cc >1 && rr>1 && cc<cols && rr<rows
            %make sure the center is the lowest graident position in 3x3
            %neighborhood - does graident nieghborhood need to be 4x4?
            subX = X(rr-1:rr+1,cc-1:cc+1,:);
            [Gmag,~,~] = imgradient3(subX);
            [~,idx] = min(Gmag(:));
            [rrSub,ccSub,~] = ind2sub(size(subX),idx);
            rrSub = rrSub-2;
            ccSub = ccSub-2;
            rr = rr+rrSub;
            cc = cc+ccSub;
            end
            
            C(1:bands+2, kk) = [squeeze(X(rr,cc,:)); cc; rr]; %could add row/col idx to feature set to keep as one 3d array %MAKE SURE CC RR ARE SUPPOSED TO FLIP HERE
            c = c+S;
            kk = kk+1;
        end
        
    end
    
    r = r+vSpacing;
end

%delete extra rows if fill pixels
if size(C,2)>kk
    C(:,kk+1:end)=[];
end

sz = size(C);
kk = sz(2);

% Now perform the clustering.  10 iterations is suggested but I suspect n
% could be as small as 2 - using 5 for speed for now
S = round(S);  % We need S to be an integer from now on

% % % %leave fill pixels clusters alone.
% % % skipK = C(1,:)==-9999;
% % % validK = 1:k;
% % % validK(skipK)=[];

%add image locations and pixel counter to the feature set
imageSize = size(X);
numRows = imageSize(1);
numCols = imageSize(2);

r = 1:numRows;
c = 1:numCols;
[c,r] = meshgrid(c,r);
pxlCounter= ones(size(r));
X = cat(3,X,r,c,pxlCounter);


for n = 1:nItr
    
    %%% ASSIGN PIXELS TO CLUSTERS
    for sp = 1:kk  % for each cluster - parfor may improve speed but require recode
        if nnz(isnan(C(1:bands,sp)))>0
            %some superpixels empty of values
            continue
        else
            % Get subimage around cluster
            rmin = max(C(bands+2,sp)-S, 1);   rmax = min(C(bands+2,sp)+S, rows);
            cmin = max(C(bands+1,sp)-S, 1);   cmax =min(C(bands+1,sp)+S, cols);
            
            subim = X(rmin:rmax, cmin:cmax, 1:bands);
            
            if numel(subim) ==  0
                warning('0 size subimage')
                continue
            else
                
                % Compute distances D between C(:,kk) and subimage
                
                D = superPixlDist(C(:, sp), subim, rmin, cmin, S, m, bands,SAMconst);
                
                
                % If any pixel distance from the cluster centre is less than its
                % previous value update its distance and label
                subd =  d(rmin:rmax, cmin:cmax);
                subl =  L(rmin:rmax, cmin:cmax);
                updateMask = D < subd;
                subd(updateMask) = D(updateMask);
                subl(updateMask) = sp;
                
                d(rmin:rmax, cmin:cmax) = subd;
                L(rmin:rmax, cmin:cmax) = subl;
            end
        end
    end
    
    %update cluster centers
    sC = regionprops(L,'centroid');
    C(:)=0;
    for j = 1:bands
        sB = regionprops(L,X(:,:,j),'MeanIntensity');
        C(j,:) = cell2mat({sB.MeanIntensity});
    end
    
    %
    % sB1 = regionprops(L,X(:,:,1),'MeanIntensity');
    % sB2 = regionprops(L,X(:,:,2),'MeanIntensity');
    % sB3 = regionprops(L,X(:,:,3),'MeanIntensity');
    % sB4 = regionprops(L,X(:,:,4),'MeanIntensity');
    % sB5 = regionprops(L,X(:,:,5),'MeanIntensity');
    % sB6 = regionprops(L,X(:,:,6),'MeanIntensity');
    % sB7 = regionprops(L,X(:,:,7),'MeanIntensity');
    % sB8 = regionprops(L,X(:,:,8),'MeanIntensity');
    % C(1:bands,:) = [sB1.MeanIntensity;sB2.MeanIntensity;sB3.MeanIntensity;sB4.MeanIntensity;...
    %     sB5.MeanIntensity;sB6.MeanIntensity;sB7.MeanIntensity;sB8.MeanIntensity];
    
    
    for sp = 1:kk
        C(bands+1,sp) = round(sC(sp).Centroid(1));
        C(bands+2,sp) = round(sC(sp).Centroid(2));
    end
    
    
    %make sure thepixels centers arn't wandering
    % % %     figure; scatter(C(10,:),C(9,:))
    % % %     title(['segmentation round ' num2str(n)])
    
    % %%%
    %     % Update cluster centres with mean values - NO NANS Allowed
    %
    % % %     validPix = l~= -1 & ~fill2D;
    % % %     idxValidPix = find(validPix);
    %     % %     imageSize = size(validPix);
    %     % %     numRows = imageSize(1);
    %     % %     numCols = imageSize(2);
    %     % %     r = 1:numRows;
    %     % %     c = 1:numCols;
    %     % %     [c,r] = meshgrid(c,r);
    %
    %     pxlIdx = 1:numRows*numCols;
    %     C(:) = 0;
    %     for ind = 1:length(pxlIdx)
    %         thisPxl = pxlIdx(ind);
    %
    %         tmp = zeros(ccData,1);
    %         rIDX = r(thisPxl);
    %         cIDX = c(thisPxl);
    %         tmp(:) = X(rIDX,cIDX,:);
    %         C(:, l(rIDX,cIDX)) = C(:, l(rIDX,cIDX)) + tmp;
    %
    %     end
    %
    %     % Divide by number of pixels in each superpixel to get mean values -
    %     % dosent count fill pixels in super pixel size
    %     for sp = 1:kk
    %         C(1:bands+2,sp) = C(1:bands+2,sp)/C(bands+3,sp);
    %     end
    %
    %     C(bands+1:bands+2,:) = round(C(bands+1:bands+2,:));
    %
    % Note the residual error, E, is not calculated because we are using a
    % fixed number of iterations
    
    fprintf(['segmentation round ' num2str(n) ' complete...\n'])
    
    
    
end
fprintf('segmentation complete!...\n')


N=length(unique(L(:)));%TEMP

%clean up the clusters
% by relabeling disjoint segments with the labels of the largest neighboring cluster.

% % % %each labeled area must be connected
% % %     %takes a long time
% % %         L = makeregionsdistinct(L);
% % %         [L, minLabel, maxLabel] = renumberregions(L);
% % %         Am = regionadjacency(L);
% % %     fprintf('cleanup complete!...\n')
% % %
% % %     N = length(Am);
% % %     toc
% % %
% changing it so only if you do clean up are these stats calculated

% N = length(unique(l(:)));

%non-contigous regions are possible, not sure if that is a problem
%though. I think adjaceny is still the same unless a fewpixels are far
%enough way to be surrounded by an region noncontigour to the main
%superpixe blob. I think my blobs are big enough that this wont be a
%problem.
%Am = regionadjacency(l); - pushing to a later step
end






function D = superPixlDist(C, im, r1, c1, S, m, bands,SAMconst)
%-- dist -------------------------------------------
%
% Usage:  D = dist(C, im, r1, c1, S, m)
%
% Arguments:   C - Cluster being considered
%             im - sub-image surrounding cluster centre
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.
%              S - grid spacing
%              m - weighting factor between colour and spatial differences.
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster centre
%
% Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% where:
% dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance
% ds = sqrt(dx^2 + dy^2)         % Spatial distance
%
% m is a weighting factor representing the nominal maximum colour distance
% expected so that one can rank colour similarity relative to distance
% similarity.  try m in the range [1-40] for L*a*b* space
%
% ?? Might be worth trying the Geometric Mean instead ??
%  Distance = sqrt(dc * ds)
% but having a factor 'm' to play with is probably handy

% Squared spatial distance
%    ds is a fixed 'image' we should be able to exploit this
%    and use a fixed meshgrid for much of the time somehow...

[rows, cols, specBands] = size(im);

%SPATIAL DISTANCE
[x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
x = x-C(bands+1);  % x and y dist from cluster centre
y = y-C(bands+2);
Dxy = sqrt(x.^2 + y.^2);

% SPECTRAL ANGLE MAPPER
clusterX =reshape(im,[rows*cols specBands]);
centerX = C(1:bands)';
DSAM = pdist2(clusterX,centerX,'cosine');
DSAM = reshape(DSAM,[rows cols]);
DSAM = DSAM*SAMconst; %get the ratio of spatial/spectral distances correct

% D = sqrt(DSAM + (dspatial2/S^2)*m^2);
D = DSAM + Dxy.*(m/S);

end