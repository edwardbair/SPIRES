function [newS] = consolidateSuperpixels(S,varargin)
% [newS] = consolidateSuperpixels(S [,tol])
%consolidate superpixel output using unique to to select a set of different
%superpixels
%
%Input
% S - output from snowSuperpixels
%Optional input
% tol - for call to uniquetol, default 0.05
%
%Output
% newS - structure S with reduced set of superpixels

p = inputParser;
addRequired(p,'S',@isstruct)
addOptional(p,'tol',.05,@(x) isscalar(x) && isnumeric(x))
parse(p,S,varargin{:});
tol = p.Results.tol;

[~,IA,~] = uniquetol(S.spectralWeightedSP,tol,'ByRows',true);
[~,IT,~] = uniquetol(S.topography,tol*3/2,'ByRows',true);
IComb = unique([IA;IT]);
newS = S;
newS.spectralWeightedSP = S.spectralWeightedSP(IComb,:);
newS.spectralSP = S.spectralSP(IComb,:);
newS.topography = S.topography(IComb,:);
newS.topoID = S.topoID;

% sort in order of descending norm
N = vecnorm(newS.spectralWeightedSP.');
[N,IA] = sort(N(:),'descend');
newS.spectralWeightedSP = newS.spectralWeightedSP(IA,:);
newS.spectralSP = newS.spectralSP(IA,:);
newS.topography = newS.topography(IA,:);
newS.normWtSP = N;

% identify likely spectrum of non-snow endmember
% criterion 1, within .05 of cosZ, so nearly flat ground
t = newS.topography(:,3)>=newS.topography(:,2)-.05 &...
    newS.topography(:,3)<=newS.topography(:,2)+.05;
% criterion 2, less than 10th percentile of those weighted reflectances
tn = newS.normWtSP<=prctile(newS.normWtSP(t),10);
possibleBackR = newS.spectralSP(tn&t,:);
newS.backReflectance = mean(possibleBackR,1,'omitnan');
newS.possibleBackR = possibleBackR;

end