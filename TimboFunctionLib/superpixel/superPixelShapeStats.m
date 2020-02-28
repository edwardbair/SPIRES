function [iReduce1,RMSE,MAE,bias,fMixedC1,medianMixC1,meadinAreaC1,iReduce2,fMixedC2,medianMixC2,meadinAreaC2] = superPixelShapeStats(L,F1,T1,F2,T2)
%SUPERPIXELSHAPESTATS Summary of this function goes here
%   L = label matrix of superpixels
%   F - fractional cover of class of interest, distributed to each pixel in
%   the superpixel.
%   T - truth mask of logical true pixels for pixels of the class of
%   interest
%
%
%
% outputs - all are fractional outputs (~range 0-1)
%ADD - statisticson the size distribution of superpixels within the image
%vs the class of interest - binary stats? for TP = SP>95% class of interest?

%HOW WELL DO SUPERPIXELS FOLLOW  BOUNDARIES?
[allSPs,idx] = unique(L);
[C1,~,~]=unique(L(T1));
[iReduce1 ]= iReduce(C1,T1);
%perfect (MAE = 0, bias = 0) would be all SPs are 100% cover of truth class
idxC = ismember(allSPs,C1);
yHat = single(F1(idx(idxC)))/100;
y = ones(size(yHat),'single');%*100;
[RMSE,MAE,bias] = rmse_mae_bias(yHat,y); %GRAPH THESE VARIABLES VS SP SIZE

A = regionprops(L,'area');
A=cell2mat({A.Area});
C1a = A(C1);
meadinAreaC1 = median(C1a);

%HOW WELL DO SUPERPIXELS SEPERATE CLASS 1 and CLASS 2?
if isempty(F2) || isempty(T2)
    iReduce2=[];
    fMixedC2=[];
    medianMixC2=[];
    meadinAreaC2=[];
    return
else
    [C2,~,~]=unique(L(T2));
    [iReduce2 ]= iReduce(C2,T2);
    C2a = A(C2);
    meadinAreaC2 = median(C2a);
    
    
    %superpixels containing class 1 and 2
    [bothC, ~, ~] = intersect(C1,C2);
    idxBoth = ismember(allSPs,bothC);
    C1frac = F1(idx(idxBoth));
    C2frac = F2(idx(idxBoth));
    
    %fraction of SP that are mixed and median mixture of classes
    [fMixedC1,medianMixC1] = fracStats(C1frac,C1);
    [fMixedC2,medianMixC2] = fracStats(C2frac,C2);
end
end

function [iRout ]= iReduce(c,t)
numSPs = length(c);
numPs = nnz(t);
iRout = (numSPs/numPs); %GRAPH VARIABLE VS SP SIZE
end

function [fMixed,medianMix] = fracStats(mixedFracs,all)
numSPs = length(all);
numMixedSPs = length(mixedFracs);
fMixed = numMixedSPs/numSPs;
medianMix = median(mixedFracs);
end
