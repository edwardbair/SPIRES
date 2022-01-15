function [newRefl] = correctAVNG(reflectance,varargin)
% [newRefl] = correctAVNG(reflectance [,threshold, default 1.6])
%correct (interpolate spatially) spurious high values of reflectance in AVIRIS-NG imagery
% reflectance should be in BIP format, use bsq2bip to change if necessary

p = inputParser;
addRequired(p,'reflectance',@isnumeric);
addOptional(p,'threshold',1.6,@(x) isnumeric(x) && isscalar(x));
parse(p,reflectance,varargin{:});
threshold = p.Results.threshold;

% indices of values above threshold
idx = find(reflectance>threshold);
% subscripts
[w,r,c] = ind2sub(size(reflectance),idx); %#ok<ASGLU>
% fill spatially
uw = unique(w).';
fprintf('%d values above threshold (%f) in %d wavelengths\n',...
    length(idx),threshold,length(uw))
newRefl = reflectance;
for k=uw
    X = squeeze(reflectance(k,:,:));
    % preserve NaNs
    t = isnan(X);
    X(t) = 0;
    mask = X>threshold;
    X = inpaintCoherent(X,mask);
    X(t) = NaN;
    newRefl(k,:,:) = X;
end
end