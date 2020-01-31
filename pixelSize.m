function [ppl,ppt,thetad] = pixelSize(R,H,p,theta_s)
% [ppl,ppt,thetad] = pixelSize(R,H,p,theta_s)
% calculate pixel sizes in along-track and cross-track directions
%
% input
%  (first 3 inputs must be in same units)
%  R - Earth radius
%  H - orbit altitude
%  p - pixel size at nadir, scalar or vector [height width]
%  theta_s - sensor zenith angle (degrees, can be vector of arguments)
%
% output (same size as vector theta_s)
%  ppl - pixel size in along-track direction
%  ppt - pixel size in cross-track direction
%  thetad - nadir sensor angle, degrees
%
% source: Dozier, J., T. H. Painter, K. Rittger, and J. E. Frew (2008),
% Time-space continuity of daily maps of fractional snow cover and albedo
% from MODIS, Advances in Water Resources, 31, 1515-1526, doi:
% 10.1016/j.advwatres.2008.08.011, equations (3)-(6).
%

% cross- and along-track pixel sizes
if isscalar(p)
    alongTrack = p;
    crossTrack = p;
else
    assert(length(p)==2,...
        'if not scalar, pixel size must be vector of length 2')
    alongTrack = p(1);
    crossTrack = p(2);
end

theta = asin(R*sind(theta_s)/(R+H));
ppl = (1/H)*(cos(theta)*(R+H)-R*sqrt(1-((R+H)/R)^2*(sin(theta)).^2));
beta = atan(crossTrack/(2*H));
ppt = (R/crossTrack)*(asin(((R+H)/R)*sin(theta+beta))-...
    asin(((R+H)/R)*sin(theta-beta))-2*beta);
thetad = radtodeg(theta);

ppl = ppl*alongTrack;
ppt = ppt*crossTrack;

end