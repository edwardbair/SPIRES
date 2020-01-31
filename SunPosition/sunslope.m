function mu = sunslope( mu0,phi0,S,A )
%sunslope cosine of illumination angle on slope
%
% mu = sunslope( mu0,phi0,S,A )
%
% Input (angles in degrees)
%   mu0, cosine of sun angle on flat surface
%   phi0, sun azimuth, degrees, +ccw from south
%   S, slope angle, degrees, from horizontal
%   A, slope azimuth, degrees, +ccw from south

mu = mu0.*cosd(S) + sqrt((1-mu0).*(1+mu0)).*sind(S).*cosd(phi0-A);
if isscalar(mu) && mu<0
    mu = 0;
else
    t = mu<0;
    mu(t) = 0;
end

end