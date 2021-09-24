function mu = sunslope( mu0,phi0,S,A )
% mu = sunslope( mu0,phi0,S,A )
%sunslope cosine of illumination angle on slope
%
% Input (angles in degrees)
%   mu0, cosine of sun angle on flat surface
%   phi0, sun azimuth, degrees, direction set by the azimuthPreference function
%   S, slope angle, degrees, from horizontal
%   A, slope azimuth, degrees, direction set by the azimuthPreference function
%
% Examples, the values for A are [-135 45] if your azimuthPreference is set
%   to convert to counter-clockwise from south, + east, - west.
%   If your azimuth preference is set to the MATLAB convention, clockwise from
%   0° north, the A should be [315 145]
%   Northern Hemisphere summer morning, Mt Blanc, NW and SE slopes
%   [declin,~,sol_lon] = EarthEphemeris(datetime('2020-07-04 10:00','TimeZone','Europe/Paris'))
%   [mu0,phi0] = sunang(45.8328,6.865,declin,sol_lon)
%   mu = sunslope(mu0,phi0,35,A)
%   Northern Hemisphere winter morning, Mt Blanc, NW and SE slopes
%   [declin,~,sol_lon] = EarthEphemeris(datetime('2020-01-04 10:00','TimeZone','Europe/Paris'))
%   [mu0,phi0] = sunang(45.8328,6.865,declin,sol_lon)
%   mu = sunslope(mu0,phi0,35,A)


[mu0,phi0,S,A] = checkSizes(mu0,phi0,S,A);
mu = mu0.*cosd(S) + sqrt((1-mu0).*(1+mu0)).*sind(S).*cosd(phi0-A);
if isscalar(mu) && mu<0
    mu = 0;
else
    t = mu<0;
    mu(t) = 0;
end

end