function f=planck(lambda,T)
% f=planck(lambda,T)
% Returns Planck radiation as a function of wavelength lambda and temperature.
%
% Wavelengths are in meters, temperatures in Kelvin.
% If the inputs are both scalars, result is a scalar.
% If one input is a vector or matrix and the other is a scalar,
% the result is a vector or matrix of the same size as the input.
% If both inputs are vectors or matrices, they must be the same size.
% If you want all combinations of 2 vectors, use meshgrid or ndgrid first.

% fundamental constants:
persistent already h k c
if isempty(already) || ~already
    C = radiationConstants;
    h = C.Planck; % Planck, J s
    k = C.Boltzmann; % Boltzmann, J/K
    c = C.SpeedOfLight; % speed of light, m/s
    already = true;
end

% if both vectors or matrices, check sizes
[lambda,T] = checkSizes(lambda,T);

assert(all(lambda(:)>=0),'wavelengths must be non-negative')
assert(all(T(:)>=0),'temperatures must be non-negative')

% Planck equation for matrices or vectors or scalars
x=h*c./(k.*lambda.*T);
f = 2*h*c^2./(lambda.^5.*expm1(x));

% error check
if any(x(:)==0) || any(isinf(x(:)))
    % check for zero or infinite wavelength, zero temperature, or negative values
    tz = isinf(x) | x==0;
    % set result for zero or infinite wavelength or zero temperature to zero
    f(tz)=0;
end
end