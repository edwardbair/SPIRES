% compares snow-covered area from 2 images, one of which is "true"
% For example, compares MODIS SCA with Landsat TM SCA.
function [bim,stats,diffimage,diffvector] = SCAcompare(C,T,scf,thresh)
% input
%   C - candidate snow-cover image
%   T - "true" snow-cover image
%       (note C & T must be same size, but T is set to NaN outside area
%       where its SCA is calculated)
%   scf - rescale factor, always < = 1 (1 is no rescaling)
%   thresh - vector of length 2, value below which C was set to zero
%       1st value is
%           tbinary - threshold for binary comparison
%       2nd value is
%           tfsca - threshold for fractional comparison
%
% output
%   bim - structure of binary images
%       tn - true negative
%       tp - true positive
%       fn - false negative
%       fp - false positive
%   stats - structure of statistics (see Wikipedia Precision_and_Recall)
%       precision
%       recall
%       specificity
%       Fstat
%       accuracy
%       mean difference
%       median difference
%       RMSE
%   diffimage - image of differences
%   diffvector - vector of differences (for pixels with snow in either image)
%
% check sizes
if ~isequal(size(C),size(T))
    error('sizes of input images must be equal')
end
%
% rescale images
if scf>1
    error('scf must be < = 1')
elseif scf<1
    Tnew = imresize(T,scf,'bilinear');
    Cnew = imresize(C,scf,'bilinear');
    if nnz(isnan(Tnew))
        Near = imresize(T,scf,'nearest');
        t = isnan(Tnew);
        Tnew(t) = Near(t);
    end
    if nnz(isnan(Cnew))
        Near = imresize(C,scf,'nearest');
        t = isnan(Cnew);
        Cnew(t) = Near(t);
    end
    T = Tnew;
    C = Cnew;
end

% apply the threshold, and also correct for possible overshoot
if isempty(thresh) || any(thresh<0) || any(isnan(thresh))
	thresh = [0 0];
elseif length(thresh)==1
    thresh(2) = 0;
end
tbinary = thresh(1);
tfsca = thresh(2);

% separate images for binary and fSCA comparisons
C(C>1) = 1;
T(T>1) = 1;
Cb = C;
Tb = T;
Cb(Cb<tbinary) = 0;
Tb(Tb<tbinary) = 0;
C(C<tfsca) = 0;
T(T<tfsca) = 0;

% compare the snow/no-snow pixels
tn = Tb==0 & Cb==0; % true negative
fp = Cb>0 & Tb==0; % false positive
fn = Cb==0 & Tb>0; % false negative
tp = Tb>0 & Cb>0; % true positive
precision = nnz(tp)/(nnz(tp)+nnz(fp));
recall = nnz(tp)/(nnz(tp)+nnz(fn));
accuracy = (nnz(tp)+nnz(tn))/(nnz(tp)+nnz(tn)+nnz(fp)+nnz(fn));
specificity = nnz(tn)/(nnz(tn)+nnz(fp));
fstat = 2*(precision*recall)/(precision+recall);

% compare the snow fractions (some are 1 or 0)
diffimage = C-T;
diffvector = diffimage(C>0 | T>0);

% binary results into structure
bim.tn = tn;
bim.tp = tp;
bim.fn = fn;
bim.fp = fp;
stats.precision = precision;
stats.recall = recall;
stats.specificity = specificity;
stats.Fstat = fstat;
stats.accuracy = accuracy;

% fractional results
diffvector(isnan(diffvector)) = [];
denom = numel(diffvector)-1;
if denom<1
    denom = 1;
end
stats.meandiff = mean(diffvector);
stats.meddiff = median(diffvector);
stats.RMSE = sqrt(sum(diffvector.^2)/denom);