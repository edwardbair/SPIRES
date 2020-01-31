function [ V ] = viewf( H, s, a )
% [ V ] = viewf( H, s, a )
%
%sky view factor for horizon circle
% (fraction of the sky open to from each point in grid)
%based on equation 7b in Dozier, J., and J. Frew (1990), Rapid calculation
%   of terrain parameters for radiation modeling from digital elevation data,
%   IEEE Trans. Geosci. Remote Sens., 28, 963-969, doi: 10.1109/36.58986.
%   (but adapted here to include Earth curvature)
%
%INPUT
%   H - sin(horizon vector) from -180 to +180 azimuths
%   s - slope, degrees
%   a - aspect, degrees, 0 south, +ccw
%
%OUTPUT
%   V - view factor for that point
%

H = double(H);
s = double(s);
a = double(a);

% azimuths of the horizon angles
phi = linspace(-pi,pi,length(H))';

% view factor
V = single(vfun(phi,asind(H),degtorad(s),degtorad(a)));

end

% integrand - eq 7b in Dozier & Frew 1990 - in that paper, H is down from
% zenith so is switched here (and H must be in radians)
function V = vfun(phi,h,slope,aspect)
h = degtorad(90-h);
cosS = cos(slope);
sinS = sin(slope);
q = cosS*sin(h).^2 + sinS.*cos(phi-aspect).*(h-sin(h).*cos(h));
V = trapz(phi,q)/(2*pi);
end