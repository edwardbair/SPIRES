function [UE] = underSegError(L,truth)
%UNDERSEGERROR Summary of this function goes here
%   Undersegmentation Error compares segment areas to measure to what
% extend superpixels flood over the ground truth segment borders. 
% For each superpixel and ground truth segment pair we count 
% the pixels in the overlap or the difference set, whichever is 
% smaller. These pixels constitute the leakage of superpixels across 
% the ground truth segments (assuming that superpixels are assigned to
%ground truth segments with highest overlap). 
% lower is better.
% A ground truth segment divides a superpixelP in aninand anoutpart. 
%oversegmentation error is the smaller error introduced by either 
%appending the out-part to the segment or by omitting the in-part of the superpixel.

%for each superpixel, calculate number of pixels inside and outside of each ground truth
%class. 
truthStats = regionprops(L,truth,'pixelvalues');
bkgrdStats = regionprops(L,~truth,'pixelvalues');


%for each ground truth object, take the minimum of each superpixels fraction (in or out) that overlaps the segment
% (both are assumed errors, assuming the Superpixel is correctly classified
spLables = unique(L(:));
spLables(spLables==0)=[];%HACK
numIN = zeros(size(spLables));
numOUT = zeros(size(spLables));
for i = 1:length(spLables)
numIN(i) = nnz(truthStats(i).PixelValues);
numOUT(i) = nnz(bkgrdStats(i).PixelValues);
end
ioCounts = [numIN numOUT];    
spE=min(ioCounts,[],2);
spError = sum(spE);
N = numel(L);
UE = spError/N;
end

