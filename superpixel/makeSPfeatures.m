function [NDSI,pixNorm,distSAM,distE] = makeSPfeatures(X,NDSIbands)
%MAKESUPERPIXELS 
%   make superpixels based on charateristics of each pixel, not just their
%   spectra.


sz = size(X);
%normalized differnce snow index
NDSI = (X(:,:,NDSIbands(1))-X(:,:,NDSIbands(2)))./(X(:,:,NDSIbands(1))+X(:,:,NDSIbands(2)));
% NDSI = normalize(NDSI(:),'range');
% NDSI = reshape(NDSI,sz(1),sz(2));

%norm (spectral vector magnitude)
p=2;
dim=3;
pixNorm = vecnorm(X,p,dim);
pixNorm = normalize(pixNorm(:),'range');
pixNorm = reshape(pixNorm,sz(1),sz(2));

%distance from "flat" reference spectra
referenceSpectra = ones(1,7).*0.5;
specX = reshape(X,sz(1)*sz(2),sz(3));

%spectral angle distance
distSAM = pdist2(referenceSpectra,specX,'cosine');
distSAM = normalize(distSAM(:),'range');
distSAM = reshape(distSAM,sz(1),sz(2));

%euclidian distance
distE = pdist2(referenceSpectra,specX,'euclidean');
distE = normalize(distE(:),'range');
distE = reshape(distE,sz(1),sz(2));
end

