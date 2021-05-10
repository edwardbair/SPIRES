function M = mieapx(cRefIn,SizPar)
% M = mieapx(cRefIn,SizPar)
% Complex angular momentum approximation to Mie scattering
% (MATLAB translation of Warren Wiscombe's FORTRAN code)
%
% Input
%   cRefIn - complex index of refraction, imaginary part can be either sign
%   SizPar - Mie size parameter (2*pi*r/wavelength)
%   (if not scalars, must be same size)
%
% Output, a structure that includes
%   omega - single-scattering albedo
%   g - asymmetry factor
%   Qext - extinction efficiency (cross-section divided by projected area)
%   Qsca - scattering efficiency
%   Qabs - absorption efficiency
%   Qpr - radiation pressure

% number of terms in CAM expansion (3 enough for 1% accuracy, 6 is the max)
nTerms = 4;
% accuracy for Qabs and Qpr, Values < 1.e-2 waste lots of computer time
Accur = 1.e-2;

% imaginary part must be positive
cRefIn = complex(real(cRefIn),abs(imag(cRefIn)));

% check sizes
[cRefIn,SizPar] = checkSizes(cRefIn,SizPar);

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

M.omega = omega;
M.g = g;
M.Qext = Qext;
M.Qabs = Qabs;
M.Qsca = Qsca;
M.Qpr = Qpr;


end