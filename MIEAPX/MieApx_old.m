function [ M, varargout ] = MieApx(radius, radiusUnits, lambda, lambdaUnits, complexIndex )
% [ M [,I] ] = MieApx(radius, radiusUnits, lambda, lambdaUnits, complexIndex )
%Mie scattering parameters for a sphere using complex angular momentum
%approximation
%
%   Input units are specified among the arguments
%   (if radius, lambda, and complexIndex are not scalars, they must be same
%   size)
% radius
% radiusUnits - any metric length ('angstrom', 'nm', 'um', 'mm', 'cm', 'm')
% lambda - wavelength of light
% lambdaUnits - any length like radiusUnits, or 'MHz', 'GHz', 'THz'
% complexIndex - of refraction, imaginary part can be positive or negative
%
% Output structure M, which includes
%   Qext - extinction efficiency
%   Qsca - scattering efficiency
%   Qabs - absorption efficiency
%   Qpr - radiation pressure efficiency
%   omega - single scattering albedo
%   g - Mie asymmetry
% Optional output structure I, which includes
%   radius - unique values only
%   lambda - unique values only
%   complexIndex - unique values only


% convert radius and lambda to meters (must be same units)
lambda = convertUnits(lambda,lambdaUnits,'m');
radius = convertUnits(radius,radiusUnits,'m');

% check and adjust input sizes if not scalars
[radius,lambda,complexIndex] = makeReplicate(radius,lambda,complexIndex);

% index of refraction of air - need to revise if we adapt this to water
% as the medium
nm = 1.000293;

% adjust wavelength for index of refraction, and calculate size parameter
waveNumber = (2*pi./lambda)*nm;
SizPar = waveNumber.*radius;
cRefIn = complexIndex/nm;

M.Qext = zeros(size(lambda));
M.Qsca = zeros(size(lambda));
M.Qabs = zeros(size(lambda));
M.Qpr = zeros(size(lambda));
M.omega = zeros(size(lambda));
M.g = zeros(size(lambda));
fn = fieldnames(M);
% CAM approximation
for idx=1:numel(cRefIn) % linear indexing even if matrix
    X = miecam(cRefIn(idx),SizPar(idx));
    for n=1:length(fn)
        M.(fn{n})(idx) = X.(fn{n});
    end
end

% optional output
if nargout>1
    I.radius = unique(radius(:));
    I.lambda = unique(lambda(:));
    I.SizPar = unique(SizPar(:));
    I.complexIndex = unique(cRefIn(:));
    varargout{1} = I;
end
end

% replicate input variables to be same size
function [R,L,N] = makeReplicate(Rin,Lin,Nin)
% which are scalars?
iscale = [isscalar(Rin) isscalar(Lin) isscalar(Nin)];
R = Rin;
L = Lin;
N = Nin;
switch nnz(iscale)
    case 0 % none is a scalar
        assert(isequal(size(R),size(L),size(N)),...
            'neither radius, lambda, nor complexIndex is scalar, they must be same size')
    case 1 % just one is a scalar
        k = find(iscale);
        switch k
            case 1 % radius is scalar
                assert(isequal(size(L),size(N)),...
                    'neither lambda nor complex index is scalar, they must be same size')
                R = Rin*ones(size(L));
            case 2 % lambda is scalar
                assert(isequal(size(R),size(N)),...
                    'neither radius nor complex index is scalar, they must be same size')
                L = Lin*ones(size(R));
            case 3 % complex index is scalar
                assert(isequal(size(R),size(L)),...
                    'neither radius nor lambda is scalar, they must be same size')
                N = Nin*ones(size(R));
        end
    case 2 % just 1 is not scalar
        k = find(~iscale);
        switch k
            case 1 % radius is not scalar
                L = Lin*ones(size(R));
                N = Nin*ones(size(R));
            case 2 % lambda is not scalar
                R = Rin*ones(size(L));
                N = Nin*ones(size(L));
            case 3 % complex index is not scalar
                R = Rin*ones(size(N));
                L = Lin*ones(size(N));
        end
    case 3 % all scalars, do nothing
end
end