function M = mieapx(cRefIn,SizPar,varargin )
% M = mieapx(cRefIn,SizPar [,useParallel] )
% Complex angular momentum approximation to Mie scattering
%
% Input
%   cRefIn - complex index of refraction, imaginary part can be either sign
%   SizPar - Mie size parameter (2*pi*r/wavelength)
%   (if not scalars, must be same size)
% Optional input
%   useParallel - true or false
%
% Output, a structure that includes
%   omega - single-scattering albedo
%   g - asymmetry factor
%   Qext - extinction efficiency (cross-section divided by projected area)
%   Qsca - scattering efficiency
%   Qabs - absorption efficiency
%   Qpr - radiation pressure

numarg = 2;
if nargin>numarg
    useParallel = varargin{1}; %#ok<NASGU>
else
    useParallel = false; %#ok<NASGU>
end

% number of terms in CAM expansion (3 enough for 1% accuracy, 6 is the max)
nTerms = 3;
% accuracy for Qabs and Qpr, Values < 1.e-2 waste lots of computer time
Accur = 1.e-1;

% check sizes
if ~isequal(size(cRefIn),size(SizPar))
    if isscalar(cRefIn)
        cRefIn = cRefIn*ones(size(SizPar));
    elseif isscalar(SizPar)
        SizPar = SizPar*ones(size(cRefIn));
    else
        error('if not scalars, size(cRefIn) must equal size(SizPar)')
    end
end

% convert matrices to vectors, and then re-convert afterward
if ~isscalar(cRefIn)
    origSize = size(cRefIn);
    cRefIn = cRefIn(:);
    SizPar = SizPar(:);
    convertSize = ~isequal(origSize,size(cRefIn));
else
    convertSize = false;
end

Qext = qexapx(cRefIn,SizPar,nTerms);
Qabs = qabapx(cRefIn,SizPar,Accur);
Qpr = qprapx(cRefIn,SizPar,Accur);
Qsca = Qext-Qabs;
omega = Qsca./Qext;
g = (Qext-Qpr)./Qsca;

% reconvert sizes
if convertSize
    Qext = reshape(Qext,origSize);
    Qabs = reshape(Qabs,origSize);
    Qpr = reshape(Qpr,origSize);
    Qsca = reshape(Qsca,origSize);
    g = reshape(g,origSize);
    omega = reshape(omega,origSize);
end

M.Qext = Qext;
M.Qabs = Qabs;
M.Qsca = Qsca;
M.Qpr = Qpr;
M.omega = omega;
M.g = g;

end