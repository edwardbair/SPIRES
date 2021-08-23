function [ mix ] = mieMixture( M, radius, density, massFraction )
% [ mix ] = mieMixture( M, radius, density, massFraction )
%Mie parameters for a mixture
%
%Input
%   M - cell vector of input Mie structures (each output from mieSphere or lookupMie)
%   radius - matrix (size(M.Qext,1+contaminants) of radii corresponding to Mie
%       structures (vectors need to be same size and all in same units, but can be any)
%   density - vector of densities of materials (same size as M)
%   massFraction - vector of mass fractions, must sum to 1.0
%
%Output
%   mix - Mie structure of the mixture
%       Qext - extinction efficiency
%       Qsca - scattering efficiency
%       Qabs - absorption efficiency
%       Qpr - radiation pressure efficiency
%       omega - single scattering albedo
%       g - Mie asymmetry parameter

% check sizes
for k=2:length(M)
    assert(isequal(size(M{1}.Qext),size(M{k}.Qext)),...
        'sizes of input Mie structure elements must be equal')
end
assert(size(radius,2)==length(M),...
    'must have one radius column for each Mie structure')
N = size(M{1}.Qext);
assert(isequal(size(M{1}.Qext,1),size(radius,1)) ||...
    size(radius,1)==1,...
    'number of elements in Mie structures and radii must be same, or there must be just a pair of radii')
assert(isequal(size(density,2),size(M,2)),'density and M must be same size')
assert(isequal(size(massFraction,2),size(M,2)),...
    'massFraction and M must have same number of columns')
if size(massFraction,1)~=1
    assert(isequal(size(massFraction,1),size(M{1}.omega,1)),...
        'massFraction and M elements must have the same number of rows')
end
% density sum within tolerance
tol = 1.e-3;
assert(max(abs(sum(massFraction,2)-1))<tol,'massFractions must sum to 1.0')
assert(all(massFraction(:)>=0 & massFraction(:)<=1),...
    'each massFraction must be in range 0-1')
S = sum(massFraction,2);
massFraction = massFraction./S;
% volume fraction
volFraction = massFraction./density;
S = sum(volFraction,2);
volFraction = volFraction./S;

% cross section follows from the ratio vol/surfA = (4/3)*pi*r^3/(4*pi*r^2)
% obviously r^2/r^3 = 1/r, but we leave in the obtuse form to allow for summing over size distributions
geomCrossSection = (3*radius.^2./radius.^3).*volFraction;

% Mie values multiplied by cross section, each row for a substance, wavelengths along rows
assert(size(geomCrossSection,2)==length(M),...
    'number of Mie structures must equal number of geometric cross sections')
gNumerator = zeros(length(M{1}.Qsca),length(M));
gDenom = zeros(size(gNumerator));
QextNumerator = zeros(size(gNumerator));
QscaNumerator = zeros(size(gNumerator));
for k=1:numel(M)
    gNumerator(:,k) = geomCrossSection(:,k).*M{k}.Qsca.*M{k}.g;
    gDenom(:,k) = geomCrossSection(:,k).*M{k}.Qsca;
    QextNumerator(:,k) = geomCrossSection(:,k).*M{k}.Qext;
    QscaNumerator(:,k) = geomCrossSection(:,k).*M{k}.Qsca;
end
gbar = sum(gNumerator,2)./sum(gDenom,2);
Qextbar = sum(QextNumerator,2)./sum(geomCrossSection,2);
Qscabar = sum(QscaNumerator,2)./sum(geomCrossSection,2);
mix.Qext = reshape(Qextbar,N);
mix.Qsca = reshape(Qscabar,N);
mix.Qabs = reshape(Qextbar-Qscabar,N);
mix.omega = reshape(Qscabar./Qextbar,N);
mix.g = reshape(gbar,N);
mix.Qpr = reshape(Qextbar-gbar.*Qscabar,N);
end