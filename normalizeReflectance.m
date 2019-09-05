function Rc=normalizeReflectance(R,topofile,solarZ,solarAzimuth)
%perform topographic correction on reflectances such that all surfaces are
%normalized to a flat and level surface
%input:
%R - reflectance cube, m x n x z - w/ m x n as height and width and z as band #
%topofile - h5 topo file from consolidateTopography
% output:
%corrected reflectance
% Reference: Tan et al 2013 http://dx.doi.org/10.1016/j.rse.2013.05.013


%b/c of Dozier's strange convention for aspect
aspect = 180 - GetTopography(topofile,'aspect'); 
slope = GetTopography(topofile,'slope');
%Tan et al eqs 1&2
c=cosd(solarZ)./(cosd(solarZ).*cosd(slope)+sind(solarZ).*...
            sind(slope).*cosd(solarAzimuth-aspect));
Rc=R.*c;