function newRefl = linearMixture(refl,fractions)
% newRefl = linearMixture(refl,fractions)
%linear mixture of spectral reflectances
%
%Input
% refl - array of N x nF reflectance values, where N is the number of
%   wavelengths and nF is the number of endmembers
% fractions - numerical vector of nF fractional values that sum to 1.0
%
%Output
% newRefl - vector of length N, reflectance of mixture

p = inputParser;
addRequired(p,'refl',@isnumeric)
addRequired(p,'fractions',@isnumeric)
parse(p,refl,fractions)
assert(size(refl,2)==length(fractions),...
    'number of columns in refl matrix must equal length of fraction vector')

% assume fractions sum to 1.0 (because fmincon has that constraint)
newRefl = zeros(size(refl,1),1);
for f=1:length(fractions)
    newRefl = newRefl+refl(:,f)*fractions(f);
end
end