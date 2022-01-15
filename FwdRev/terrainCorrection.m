function [newRefl] = terrainCorrection(Prescription,directRefl)
% [newRefl] = terrainCorrection(Prescription,directRefl)
%correct spectral reflectance for terrain
%Input
% Prescription - original prescription from which direct reflectance was
%   calculated
% directRefl - direct reflectance at angle Prescription.Illumination.muS
%
%Output
% newRefl - vector of reflectance values accounting for direct and diffuse
%   illumination over terrain

p = inputParser;
addRequired(p,'Prescription',@isstruct);
addRequired(p,'directRefl',@isnumeric);
parse(p,Prescription,directRefl);

Illumination = Prescription.Illumination;
% check Illumination values
assert(~isempty(Illumination.cosZ) && ~isempty(Illumination.muS),...
    'Illumination structure must not have empty cosZ or muS')
assert(Illumination.cosZ>0 && Illumination.cosZ<=1,...
    'Illumination structure cosZ  must be >0 and <=1')
assert(Illumination.muS>=0 && Illumination.muS<=1,...
    'Illumination structure muS  must be >=0 and <=1')
assert(Illumination.viewF>0 && Illumination.viewF<=1,...
    'Illumination view facture must be >0 and <=1')

% calculate diffuse fraction
diffuseFraction = getAtmosProp(Prescription.RadiativeTransfer.diffF,...
    Prescription.Spectrum.wavelength,Illumination.elevation,Illumination.cosZ);

% calculate diffuse reflectance with same prescription but angles
% eliminated
tmpP = Prescription;
tmpP.RadiativeTransfer.calcTerrain = false;
tmpP.cosZ = [];
tmpP.muS = [];
diffuseRefl = SPIReS_fwd(tmpP);
newRefl = directRefl.*(1-diffuseFraction)*Illumination.muS/Illumination.cosZ +...
    diffuseRefl.*diffuseFraction*Illumination.viewF;
end