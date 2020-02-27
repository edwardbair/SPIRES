function [L] = SAMSuperPixels(superPixelSize,X,m,SAMconst,nItr,validPixels)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% L- label matrix of the scene
%       positive lables - possible clouds
%       0 label - not cloud
%       -1 - background/fill

% SPECTRAL ANGLE MAPPER SUPER PIXELS
[nR, nC,~] = size(X);
k = round((nR.*nC)./superPixelSize); %number of desired superPixels
[L,~,~] = slicMultiSpec(X, k, m, nItr,SAMconst);
L(~validPixels)=0;
end

