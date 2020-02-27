function  [featureSet] = calcSuperPixelFeatures(L,varargin)
%[rmCFMask] = CFMaskRemove(L,GFBTO,w653,snow,notSnow)
%pixels to remove from CFMask using gabor filterbank output and superpixels
%   Detailed explanation goes here
% L- labeled regions of the superpixels - positive integers. 0 - background
% pixels with no charateristics calculated.
%
% varargin - pixel values for feature of interest that fractional coverage
% of superpixel will be calculated for
%
%output
%   featureSet - superpixel charateristcs distributed to each pixel in the
%   superpixel as percentages in uint8 integers (0-100)


assert(exist('varargin','var')==1,'must include at least one feature for superpixel calculations')

for i = 1:length(varargin)
    assert(islogical(varargin{i}),'pixel features must be logical matrices')
end

checkSuperPixels = unique(L(~L==0));
boxes=regionprops(L,'BoundingBox');

featureSet = zeros([size(L) length(varargin)],'uint8');

for sp = 1:length(checkSuperPixels)  % for each cluster - parfor may improve speed but require recode
    
    thisSP = checkSuperPixels(sp);
    
    % Get subimage around cluster [left, top, width, height]
    BB = boxes(thisSP);
    rmin = ceil(BB.BoundingBox(2));
    rmax = rmin+BB.BoundingBox(4)-1;
    cmin = ceil(BB.BoundingBox(1));
    cmax = cmin + BB.BoundingBox(3)-1;
    subL = L(rmin:rmax,cmin:cmax);
    if numel(subL) ==  0
        warning('0 size subimage')
        continue
    end
    subSP = subL == thisSP;
    [row,col]=find(subSP);
    numPxls = nnz(subSP);
    
    %calulate superpixel fractional cover of true
    
    for k = 1:length(varargin)
        subFeatureSet= featureSet(rmin:rmax,cmin:cmax,k);
        F = varargin{k};
        subF = F(rmin:rmax,cmin:cmax);
        
        %round up to catch small inclusions of pixels less than 1% of SP
        subFeatureSet(subSP)=ceil((nnz(subF&subSP)/numPxls)*100);
        
        %fill pixels in full mask with superpixel fractional cover of true pixels
        featureSet(rmin:rmax,cmin:cmax,k) = subFeatureSet;
    end
    
    
end
end

