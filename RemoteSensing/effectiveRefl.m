function R = effectiveRefl(R0,Rdif,mu0,muS,diffuseFraction)
%effective reflectance as a function of illumination angles
%
%Input
%   R0 - reflectance (Directional-Hemispherical) of surface at illumination
%       angle on slope
%   Rdif - diffuse reflectance (Hemispherical-Hemispherical) of surface
%       (equivalent to R0 integrated over all illumination angles)
%   mu0 - cosine of illumination angle on flat surface
%   muS - cosine of illumination angle on slope
%   diffuseFraction - fraction (0-1) of incoming radiation that is diffuse
%
%Output
%   R - effective reflectance, as seen from above without knowing the slope
%       or illumination geometry

[R0,Rdif,mu0,muS,f] = checkSizes(R0,Rdif,mu0,muS,diffuseFraction);
R = (1-f).*R0.*muS./mu0 + f.*Rdif;
end