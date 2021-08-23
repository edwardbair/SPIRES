function [ M ] = mieSphere(radius, radiusUnits, lambda, lambdaUnits, varargin )
% [ M ] = mieSphere(radius, radiusUnits, lambda, lambdaUnits, ... )
%Basic Mie scattering parameters for a sphere, calculated either with
%Wiscombe's complex angular momentum approximation or the translation of his
%MIEV0 FORTRAN code, depending on size parameter x = 2 pi r / wavelength
%
%   Input units are specified among the arguments
% radius
% radiusUnits - any metric length ('angstrom', 'nm', 'um' or 'mum', 'mm', 'cm', 'm')
% lambda - wavelength of light
% lambdaUnits - any length like radiusUnits, or 'MHz', 'GHz', 'THz'
%
% Optional inputs, name-value pairs
%   Either 'substance' or 'refindex' must be specified
%   'substance' - either 'ice', 'snow', 'water', 'dust', or 'soot' to calculate
%       refractive index
%       or
%   'refindex' - complex number(s), imaginary part either positive or negative,
%       can be a vector or matrix the same size as lambda, or a constant
%   'mediumRI' - index of refraction in medium, default air 1.000293,
%       either scalar or same length as refindex
%   'CAMthreshold' - scalar to specify critical Mie parameter below which
%       the regular Mie approach and above which the complex angular
%       momentum method is used, default 20, could go as high as 100
%
% Output structure M, which includes
%   Qext - extinction efficiency
%   Qsca - scattering efficiency
%   Qabs - absorption efficiency
%   Qpr - radiation pressure efficiency
%   omega - single scattering albedo
%   g - Mie asymmetry parameter

% defaults
defaultSubstance = '';
defaultRefIndex = [];
% parse inputs
p = inputParser;
positiveValidation = @(x) isnumeric(x) && all(x(:)>0);
addRequired(p,'radius',positiveValidation)
addRequired(p,'radiusUnits',@ischar)
addRequired(p,'lambda',positiveValidation)
addRequired(p,'lambdaUnits',@ischar)
addParameter(p,validatestring('substance',{'subs','substance'}),...
    defaultSubstance,@ischar)
addParameter(p,validatestring('refindex',{'refi','refindex'}),...
    defaultRefIndex,@isnumeric)
addParameter(p,validatestring('mediumri',{'med','mediumr','mediumri'}),...
    1.000293,@(x) isnumeric(x))
addParameter(p,validatestring('camthreshold',{'camt','camthresh','camthreshold'}),...
    20,@(x) isnumeric(x) && isscalar(x))
parse(p,radius,radiusUnits,lambda,lambdaUnits,varargin{:})
radius = p.Results.radius;
radiusUnits = p.Results.radiusUnits;
lambda = p.Results.lambda;
lambdaUnits = p.Results.lambdaUnits;
% either substance or refractive index must be specified
assert(xor(isempty(p.Results.substance),isempty(p.Results.refindex)),...
    '''substance'' or ''refindex'' must be specified, but not both')

% convert radius and lambda to meters (needed to calculate size parameter)
lambda = convertLengthUnits(lambda,lambdaUnits,'m');
radius = convertLengthUnits(radius,radiusUnits,'m');

% index of refraction, either specified or by substance
if ~isempty(p.Results.substance)
    complexIndex = refractiveIndex(lambda,p.Results.substance,'m');
else
    if isreal(p.Results.refindex)
        complexIndex = complex(p.Results.refindex,0);
    else
        complexIndex = p.Results.refindex;
    end
end
[complexIndex,lambda,radius,nm] =...
    checkSizes(complexIndex,lambda,radius,p.Results.mediumri);

% size parameter and refractive index in medium with index nm
waveNo = (2*pi./lambda).*nm;   % wavenumber in medium nm
xx = waveNo.*radius;        % size parameter
N = complexIndex./nm;        % relative refractive index

% simple if both scalars
if isscalar(xx) && isscalar(N)
    if xx>=p.Results.camthreshold && real(N)>1 %#ok<BDSCI>
        M = mieapx(N,xx);
    else
        M = mieQvals(N,xx);
    end
    return
end

% which values can use mieapx?
sizeProblem = size(xx);
okay = xx>=p.Results.camthreshold & real(N)>1;
if all(okay)
    M = mieapx(N,xx);
elseif all(~okay)
    M = mieQvals(N,xx);
else
    Ma = mieapx(N(okay),xx(okay));
    Mm = mieQvals(N(~okay),xx(~okay));
    % put the two structure together
    fn = fieldnames(Ma);
    for f=1:length(fn)
        M.(fn{f}) = zeros(sizeProblem);
        M.(fn{f})(okay) = Ma.(fn{f});
        M.(fn{f})(~okay) = Mm.(fn{f});
    end
end
end