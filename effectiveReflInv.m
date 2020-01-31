function R0 = effectiveReflInv(Reff,mu0,muS,secS)
%inverse of effectiveRefl, estimate surface reflectance from effective
%reflectance
%
%Input
%   Reff - apparent reflectance (Directional-Hemispherical) as seen from
%       above without knowing slope or illumination geometry
%   mu0 - cosine of illumination angle on flat surface
%   muS - cosine of illumination angle on slope
%   secS - secant of slope angle (i.e. 1/cos)
%
%Output
%   R0 - reflectance on slope at illumination angle acos(muS)

R0 = Reff.*mu0.*secS./muS;
end

