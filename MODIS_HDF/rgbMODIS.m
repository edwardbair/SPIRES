function [ rgb ] = rgbMODIS( I,bands,varargin )
% [ rgb ] = rgbMODIS( I,bands )
%create RGB image from MODIS 7-band reflectance image
%
% Input
%   I - MODIS image, all 7 bands, for example from GetMOD09GA
%   bands - 3-element vector of bands to use
%       for example SWIR NIR VIS = [6 2 4]
%           traditional false color NIR Red Green [2 1 4]
%           true color Red Green Blue [1 4 3]
% Optional input
%   'noeq' - do not equalize histograms (default is to run histeq)

X = double(cat(3,I(:,:,bands(1)),I(:,:,bands(2)),I(:,:,bands(3))));
X(isnan(X)) = 0;
if ~(nargin>2 && strcmp(varargin{1},'noeq'))
    for m=1:size(X,3)
        X(:,:,m) = histeq(X(:,:,m));
    end
end
X(X>1) = 1;
X(X<0) = 0;
rgb = uint8(X*double(intmax('uint8')));
end