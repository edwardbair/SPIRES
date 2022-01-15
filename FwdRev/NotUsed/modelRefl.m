function [refl] = modelRefl(prescription)
%assemble modeled reflectance based on snow properties, and possibly terrain
%and diffuse fraction
%
%Input
% prescription

p = inputParser;
addRequired(p,'prescription',@isstruct)
parse(p,prescription);
if isempty(p.Results.diffuseFraction)
    refl = SPIReS_fwd(prescription);
else
    dF = p.Results.diffuseFraction;
    R1 = SPIReS_fwd(prescription);
    cosineRatio = prescription.Illumination.muS/prescription.Illumination.cosZ;
    prescription.Illumination.cosZ = [];
    prescription.Illumination.muS = [];
    R2 = SPIReS_fwd(prescription);
    refl = R1.*(1-dF)*cosineRatio+R2.*dF*prescription.Illumination.viewF;
end