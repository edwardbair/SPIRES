function [C] = calcSuperPixelMeans(L,X)
%CALCSUPERPIXELMEANS 
%calculate mean feature values in the superpixels
%
%input
%   L - MxN matrix of labels, positive integers for which superpixel each pixel belongs
%   to, 0 values are background/fill pixels that are skipped
%   X - MxNx(number of features) for which means shoudl be calculated
% output
%   C - table of mean values for each superpixel, ech superpixel is a row
%   row number is same as label number in L.

numBands = size(X,3);
k = length(unique(L(L~=0)));
C = zeros(k,numBands);

for j = 1:numBands
    sB = regionprops(L,X(:,:,j),'MeanIntensity');
    C(:,j) = cell2mat({sB.MeanIntensity});
end

end

