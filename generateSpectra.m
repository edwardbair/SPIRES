function Rmix = generateSpectra(R,backgroundR,fSCA, varargin)
% Rmix = generateSpectra(R,backgroundR,fSCA, varargin)
%generate mixed-pixel spectra from calculation of 100% snow spectrum
%
%Input
%   R - input spectrum, calculated for pixel with fSCA=1
%   backgroundR - reflectance of pixel without snow in same bands as R
%   fSCA - fractional snow cover (0-1)
%Optional input, name-value pairs
%   'noise' - random error (fraction) in mixture (std dev of normally distributed
%       error)
%   'bias' - bias (fraction) in mixture, between +/- 0.2;
%
%Output
%   Rmix - reflectance of mixed pixel, with or without random error

p = inputParser;
addRequired(p,'R',@isnumeric)
addRequired(p,'backgroundR',@isnumeric)
addRequired(p,'fSCA',@isnumeric)
addParameter(p,'noise',0,@isnumeric)
addParameter(p,'bias',0,@isnumeric);
parse(p,R,backgroundR,fSCA,varargin{:})

% size check
[R,bR] = checkSizes(p.Results.R,p.Results.backgroundR);

% mixed pixel
Rmix = p.Results.fSCA*R+(1-fSCA)*bR;
if any(Rmix(:)<0 | Rmix(:)>1)
    warning('some Rmix outside range 0-1, reset')
    Rmix(Rmix<0) = 0;
    Rmix(Rmix>1) = 0.99;
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
assert(abs(p.Results.bias)<=0.2)
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