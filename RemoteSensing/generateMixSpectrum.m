function Rmix = generateMixSpectrum(R,fmix, varargin)
% Rmix = generateMixSpectrum(R,fmix, varargin)
%generate mixed-pixel spectrum from matrix of spectra of endmembers
%
%Input
%   R - input spectra matrix of size [nBands, nMembers]
%   fmix - fractional coverage of endmembers (automatically adjusted to sum
%       to 1.0)
%Optional input, name-value pairs
%   'noise' - random error (fraction) in mixture (std dev of normally distributed
%       error)
%   'bias' - bias (fraction) in mixture, between +/- 0.2;
%
%Output
%   Rmix - reflectance of mixed pixel, with or without random error

p = inputParser;
RValidation = @(x) isnumeric(x) && ismatrix(x) &&...
    all(x(:)>=0) && all(x(:)<=1);
MValidation = @(x) isnumeric(x) && (isrow(x) || iscolumn(x)) &&...
    all(x(:)>=0) && all(x(:)<=1);
NValidation = @(x) isnumeric(x) && isscalar(x);
addRequired(p,'R',RValidation)
addRequired(p,'fmix',MValidation)
addParameter(p,'noise',0,NValidation)
addParameter(p,'bias',0,NValidation)
parse(p,R,fmix,varargin{:})

% size check
assert(isequal(size(p.Results.R,2),length(p.Results.fmix)),...
    'size(R) = [%g %g], length(fmix)=%g should be same as #columns in R',...
    size(p.Results.R,1),size(p.Results.R,2),length(p.Results.fmix))

% set sum(fmix) to 1.0
fmix = p.Results.fmix/sum(p.Results.fmix);

% mixed pixel
R = p.Results.R;
Rmix = zeros(length(R),1);
for k=1:length(fmix)
    Rmix = Rmix+R(:,k)*fmix(k);
end

% add random error
assert(p.Results.noise>=0 && p.Results.noise<1,...
    'randErr must be >=0 and <1')
if p.Results.noise~=0
    rng('shuffle')
    errValue = (randn(size(Rmix)))*p.Results.noise;
    newRmix = Rmix.*(1-errValue);
else
    newRmix = Rmix;
end

% add bias
assert(abs(p.Results.bias)<=0.2,'bias must be in range +/- 0.2')
if p.Results.bias~=0
    newRmix = newRmix*(1+p.Results.bias);
end

% make sure we stay within 0-1
t0 = newRmix<0;
t1 = newRmix>=1;
if nnz(t0)
    newRmix(t0) = Rmix(t0);
end
if nnz(t1)
    newRmix(t1) = Rmix(t1);
end
Rmix = newRmix;
end