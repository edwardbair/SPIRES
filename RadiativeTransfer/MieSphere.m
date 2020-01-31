function [ M ] = MieSphere(radius, radiusUnits, lambda, lambdaUnits, varargin )
% [ M ] = MieSphere(radius, radiusUnits, lambda, lambdaUnits, ... )
%Basic Mie scattering parameters for a sphere, calculated either with
%MatScat or Wiscombe's complex angular momentum approximation, depending on
%size
%
%   Input units are specified among the arguments
% radius
% radiusUnits - any metric length ('angstrom', 'nm', 'um' or 'mum', 'mm', 'cm', 'm')
% lambda - wavelength of light
% lambdaUnits - any length like radiusUnits, or 'MHz', 'GHz', 'THz'
%
% Optional inputs, name-value pairs
%   'useParallel' - true to use a parfor loop for the calculations, default
%       false
%   Either 'substance' or 'refindex' must be specified
%   'substance' - either 'ice', 'snow', 'water', 'dust', or 'soot' to calculate
%       refractive index
%       or
%   'refindex' - complex number(s), imaginary part either positive or negative,
%       can be a vector or matrix the same size as lambda, or a constant
%
% Output structure M, which includes
%   Qext - extinction efficiency
%   Qsca - scattering efficiency
%   Qabs - absorption efficiency
%   Qpr - radiation pressure efficiency
%   omega - single scattering albedo
%   g - Mie asymmetry parameter

% defaults
defaultParallel = false;
defaultSubstance = '';
defaultRefIndex = [];
% parse inputs
p = inputParser;
positiveValidation = @(x) isnumeric(x) && all(x(:)>0);
addRequired(p,'radius',positiveValidation)
addRequired(p,'radiusUnits',@ischar)
addRequired(p,'lambda',positiveValidation)
addRequired(p,'lambdaUnits',@ischar)
addParameter(p,'useparallel',defaultParallel,@islogical)
addParameter(p,'substance',defaultSubstance,@ischar)
addParameter(p,'refindex',defaultRefIndex,@isnumeric)
parse(p,radius,radiusUnits,lambda,lambdaUnits,varargin{:})
radius = p.Results.radius;
radiusUnits = p.Results.radiusUnits;
lambda = p.Results.lambda;
lambdaUnits = p.Results.lambdaUnits;
useParallel = p.Results.useparallel;
% either substance or refractive index must be specified
assert(xor(isempty(p.Results.substance),isempty(p.Results.refindex)),...
    '''substance'' or ''refindex'' must be specified, but not both')

% convert radius and lambda to meters (needed to calculate size parameter)
lambda = convertUnits(lambda,lambdaUnits,'m');
radius = convertUnits(radius,radiusUnits,'m');

% index of refraction, either specified or by substance
if ~isempty(p.Results.substance)
    complexIndex = RefractiveIndex(lambda,p.Results.substance,'m');
else
    if isreal(p.Results.refindex)
        complexIndex = complex(p.Results.refindex,0);
    else
        complexIndex = p.Results.refindex;
    end
    % imaginary part positive
    t = imag(complexIndex)<0;
    if nnz(t)
        complexIndex(t) = conj(complexIndex(t));
    end
end
[complexIndex,lambda,radius] = checkSizes(complexIndex,lambda,radius);

% index of refraction of air - need to revise if we adapt this to water
% as the medium
nm = 1.000293;

% size parameter and refractive index in medium with index nm
waveNo = 2*pi./lambda*nm;   % wavenumber in medium nm
xx = waveNo.*radius;        % size parameter
N = complexIndex/nm;        % relative refractive index
% this routine doesn't deal with real(N)<1
t = real(N)<1;
if nnz(t)
    warning('%d values have real(N)<1 (min %f), which this routine doesn''t deal with so adjusted to %f',...
        nnz(t),min(real(N)),1+10*eps);
    N(t) = complex(1+10*eps,imag(N(t)));
end

% simple if both scalars
if isscalar(xx) && isscalar(N)
    if xx>10 && real(N)>1 %#ok<BDSCI>
        M = mieapx(N,xx);
    else
        M = MieStuff(radius,lambda,N,nm);
    end
    return
end

% which values can use mieapx?
sizeProblem = size(xx);
okay = xx>10 & real(N)>1;
if all(okay)
    M = mieapx(N,xx);
elseif all(~okay)
    M = MieStuff(radius,lambda,N,nm,useParallel);
else
    Ma = mieapx(N(okay),xx(okay));
    Mm = MieStuff(radius(~okay),lambda(~okay),N(~okay),nm,useParallel);
    % put the two structure together
    fn = fieldnames(Ma);
    for f=1:length(fn)
        M.(fn{f}) = zeros(sizeProblem);
        M.(fn{f})(okay) = Ma.(fn{f});
        M.(fn{f})(~okay) = Mm.(fn{f});
    end
end
end