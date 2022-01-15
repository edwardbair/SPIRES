function T=planckinv(lambda,L)
% T=planckinv(lambda,L)
% Returns Planck blackbody temperature as a function of wavelength and radiance.
% Wavelength lambda in meters, radiance L in W m^-3 sr^-1 (same as output
%   from planck.m)
% If the inputs are both scalars, result is a scalar.
% If one input is a vector or matrix and the other is a scalar,
% the result is a vector or matrix of the same size as the input.
% If both inputs are vectors or matrices, they must be the same size.
% If you want all combinations of 2 vectors, use meshgrid or ndgrid first.

% fundamental constants:
persistent already h k c
if isempty(already) || ~already
    C = radiationConstants;
    h = C.Planck; % Planck, J/s
    k = C.Boltzmann; % Boltzmann, J/K
    c = C.SpeedOfLight; % speed of light, m/s
    already = true;
end

[lambda,L] = checkSizes(lambda,L);
assert(all(L(:)>0),'L must be positive')
assert(all(lambda(:)>0),'lambda must be positive')
assert(~any(isinf(L(:))),'L cannot be infinite');
assert(~any(isinf(lambda(:))),'lambda cannot be infinite');

y = 2*h*c*c./(lambda.^5.*L);
x = log1p(y);
T = h*c./(k.*lambda.*x);
end