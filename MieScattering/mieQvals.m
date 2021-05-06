function M = mieQvals(cRefIn,SizPar)
% M = mieQvals(cRefIn,SizPar)
%This Mie scattering program complements mieapx in providing just the Mie
%parameters needed to calculate fluxes, i.e. extinction, scatter,
%absorption, and radiation pressure efficiencies, along with single
%scattering albedo and Mie asymmetry parameter.
%It's a MATLAB translation of Warren Wiscombe's miev0 code, also taking
%advantage of Scott Prahl's C and Python codes.
%
%Inputs can be scalars, vectors, matrices, or N-dimensional objects, but
%any non-scalar entries must be the same size
% cRefIn - complex index of refraction with imaginary part negative (routine
%   checks and fixes this)
% SizPar - Mie size parameter
%
% Output structure M, which includes
%   Qext - extinction efficiency
%   Qsca - scattering efficiency
%   Qabs - absorption efficiency
%   Qpr - radiation pressure efficiency
%   omega - single scattering albedo
%   g - Mie asymmetry parameter
%
%For size parameter x >= ~20, mieapx is recommended instead

% imaginary part must be negative
if any(imag(cRefIn(:))~=0)
    cRefIn = complex(real(cRefIn),-abs(imag(cRefIn)));
end

[m,x] = checkSizes(cRefIn,SizPar);
% convert matrices to vectors, and then re-convert afterward
if ~isscalar(m)
    origSize = size(m);
    x = x(:);
    m = m(:);
    convertSize = ~isequal(origSize,size(m));
else
    [qext,qsca,~,g] = mieScalar(m,x); % qback not used
    M.Qext = qext;
    M.Qsca = qsca;
    M.Qabs = qext-qsca;
    M.omega = qsca/qext;
    M.g = g;
    M.Qpr = qext-g*qsca;
    return
end

% general case, m & x are vectors of the same size
qext = zeros(size(x));
qsca = zeros(size(x));
g = zeros(size(x));
for k=1:length(x)
    % mieScalar calculates qback, but it's not used (but left in the
    % mieScalar code just in case)
    [qext(k),qsca(k),~,g(k)] = mieScalar(m(k),x(k));
end
if convertSize
    qext = reshape(qext,origSize);
    qsca = reshape(qsca,origSize);
    g = reshape(g,origSize);
end
M.Qext = qext;
M.Qsca = qsca;
M.Qabs = qext-qsca;
M.omega = qsca./qext;
M.g = g;
M.Qpr = qext-g.*qsca;
end