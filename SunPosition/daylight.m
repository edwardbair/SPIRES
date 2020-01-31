function omega = daylight(doy, latitude)
% omega = daylight(doy, latitude)
%
% length of day, does not consider refraction
%
% Input
%   doy - calendar day of year (1-366)
%   latitude - in degrees, + in northern hemisphere
%   (if both arguments are not scalars, they need to be same size)
%
% Output
%   omega - length of daylight in degrees of rotation (i.e. 360 = 24 hrs)

if ~isscalar(doy) && ~isscalar(latitude)
    assert(isequal(size(doy),size(latitude)),...
        'if not scalars, arguments must be of same size')
end

declination = declination_simple(doy);
arg = -tand(latitude).*tand(declination);
% reset argument for total darkness (>1) or total sunlight(<-1)
t = abs(arg)>1;
if nnz(t)
    arg(t) = sign(arg(t));
end
omega = 2*acosd(arg);

end