function [wt] = getAtmosProp(F,wavelength,elevation,cosZ)
% [wt] = getAtmosProp(F,wavelength,elevation,cosZ)
%   get atmospheric property, either transmittance weight or diffuse fraction
%
% Input
%   F - either Ftrans or Fdf (from buildAtmosWeightF)
%   wavelength - vector or matrix, same units for wavelength as in F
%   elevation - meters, scalar or same size as wavelength
%   cosZ - cosine of illumination angle, scalar or same size as wavelength

[w,z,mu0] = checkSizes(double(wavelength),double(elevation),double(cosZ));
wt = F(w,z,mu0);
end

