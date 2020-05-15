function [Rc,c]=normalizeReflectance(R,slope,aspect,solarZ,solarAzimuth)
%perform topographic correction on reflectances such that all surfaces are
%normalized to a flat and level surface
%input:
%R - reflectance cube, m x n x z - w/ m x n as height and width and z as band #
%slope and aspect, e.g. from consolidateTopography
% output:
%corrected reflectance
% Reference: Tan et al 2013 http://dx.doi.org/10.1016/j.rse.2013.05.013
% c - correction factor

%b/c of Dozier's strange convention for aspect
aspect = 180 - aspect;
%slope = GetTopography(topofile,'slope');
%Tan et al eqs 1&2
ic=(cosd(solarZ).*cosd(slope)+sind(solarZ).*...
            sind(slope).*cosd(solarAzimuth-aspect));

c=cosd(solarZ)./ic;
Rc=R.*c;

%fix negative values
Rc(Rc<0)=0;
% fix values > 1
% Rc(Rc>1)=1;
