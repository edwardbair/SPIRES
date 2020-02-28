function [recall,classL] = boundaryRecall(L,d,truth,classFlag)
%BOUNDARYRECALL Summary of this function goes here
%   fraction of ground truth edges that fall within a certain distance of a least one superpixel boundary.
% L = superpixelLabels
% d = distance from boundary that counts as a true positive boundary
% truth - image with truth class (function calculates the boudary
% internally
%higher is better
% true positive boundary pixels in S with respect to ground truth.


%flag to limit to only the class of interest for if it dosent cover the
%whole scene?

%assume superpixels w/majority pixel class are correctly classifiable as
%that class, how well does the exterior boundary of the class mapped by the
%superpixel?

if classFlag
    allSP = unique(L);
    validSP = unique(L(truth));
    %lables in L but not valid becuase no truth class cover
    removeSP = setdiff(allSP,validSP);
    removeSP(removeSP==0)=[];
    
    s = regionprops(L,'PixelIdxList');
    for k = 1:length(removeSP)
        rm = removeSP(k);
        idx = s(rm).PixelIdxList;
        L(idx) = 0;
    end
    
    %labels in L but not valid becuase less than 50% truth class cover
    p = regionprops(L,truth,'pixelvalues');
    s2 = regionprops(L,'PixelIdxList');
    spLables = unique(L);
    spLables(spLables==0)=[];
    fracCover= zeros(size(spLables));
    for i = 1:length(spLables)
        fracCover(i) = nnz(p(i).PixelValues)/length(p(i).PixelValues);
        if fracCover(i)<1 && spLables(i)~=0
            rm = spLables(i);
            idx = s2(rm).PixelIdxList;
            L(idx) = 0;
        end
    end   
    %merge all superpixels into a single class superpixel
    L = L~=0;   
end


%4-way connected boundary for superpixels
spBM= boundarymask(L,4);

%enlarge truth boundaries by d pixels as close enough to count as a TP
%boundary (also 8-way connected)
truthBM = boundarymask(truth,8);
SE = strel('square',d);
truthBM = imdilate(truthBM,SE);

%True Positives (TP)Number of boundary pixels in G for whose exista boundary pixel in B in ranged.
TP = spBM & truthBM;
%False Negatives (FN)Number of boundary pixels in G for whose doesnot exist a boundary pixel in B in ranged.
FN = spBM & ~truthBM;

% Boundary Recall R = TP/(TP+FN)
recall = nnz(TP)/(nnz(TP)+nnz(FN));
classL=L;
end

