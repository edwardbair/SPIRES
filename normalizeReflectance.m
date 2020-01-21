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
ic=(cosd(solarZ).*cosd(slope)+sind(solarZ).*...
            sind(slope).*cosd(solarAzimuth-aspect));

c=cosd(solarZ)./ic;
Rc=R.*c;
%limit corrections


%fix negative values
Rc(Rc<0.001)=0.001;
%fix values > 1
Rc(Rc>=1)=0.999;