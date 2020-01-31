function R = effectiveRefl(R0,mu0,muS,cosS)
%effective reflectance as a function of slope and illumination angles
%
%Input
%   R0 - reflectance (Directional-Hemispherical) of surface at illumination
%       angle on slope
%   mu0 - cosine of illumination angle on flat surface
%   muS - cosine of illumination angle on slope
%   cosS - cosine of slope angle
%
%Output
%   R - effective reflectance, as seen from above without knowing the slope
%       or illumination geometry

R = R0.*cosS.*muS./mu0;
end

