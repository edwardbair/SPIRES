function omega = sunhorz( latitude, longitude, declination, varargin )
% omega = sunhorz( lat, lon, declin, [P,T] )
%sunhorz longitudes of sun at sunrise & sunset with option to account for
%       refraction
%  Input, all in degrees, all scalar
%   latitude
%   longitude
%   declination
%
%  Variable scalar input arguments, must be 2 to account for refraction
%   P - atmospheric pressure, kPa
%   T - temperature, Kelvin
%
%  Output:
%   omega, vector of length 2, solar longitudes at sunrise & sunset
%   (east of location at sunrise, west at sunset)

global LAT DECLIN PRESS TEMP ZEROFLAG

assert(isscalar(latitude) && isscalar(longitude) && isscalar(declination),...
    'all inputs must be scalars')
arg = -tand(latitude).*tand(declination); % now scalars, but may incorporate vectors later
% reset argument for total darkness (>1) or total sunlight(<-1)
if abs(arg)>1
    arg = sign(arg);
end

optarg = size(varargin,2);
assert(optarg == 0 || optarg ==2, 'must have 0 or 2 optional arguments')
if optarg == 0
    % standard method, ignoring atmospheric refraction
    o = acosd(arg);
else
    % numerical method accounting for atmospheric refraction
    PRESS = varargin{1};
    TEMP = varargin{2};
    ZEROFLAG = false;
    LAT = latitude;
    DECLIN = declination;  
    o = fzero(@psunangr,acosd(arg));
end

omega = [longitude+o longitude-o];
t = omega>180;
omega(t) = omega(t)-180;
t = omega<-180;
omega(t) = abs(180+omega(t));

end