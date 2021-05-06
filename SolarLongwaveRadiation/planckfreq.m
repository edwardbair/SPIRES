function f=planckfreq(nu,T)
% f=planckfreq(nu,T)
%
% Returns Planck radiation as a function of frequency and temperature.
% Frequency nu in Hertz, temperature T in Kelvin.
% If either frequency or temperature is a vector, the return values are a
% vector.
% If both frequency and temperature are vectors or matrices, they must be
% the same size.
% If you want all combinations of two vectors, use meshgrid first.

% fundamental constants:
persistent already h k c
if isempty(already) || ~already
    C = radiationConstants;
    h = C.Planck; % Planck, J/s
    k = C.Boltzmann; % Boltzmann, J/K
    c = C.SpeedOfLight; % speed of light, m/s
    already = true;
end

[nu,T] = checkSizes(nu,T);

assert(all(nu(:)>=0),'all frequencies must be non-negative')
assert(all(T(:)>=0),'all temperatures must be non-negative')

x=h*nu./(k*T);
f=2*h*nu.^3./(c^2.*expm1(x));

if any(isinf(x(:))) || any(x(:)==0)
    tz = isinf(x) | x==0;
    f(tz) = 0;
end

end