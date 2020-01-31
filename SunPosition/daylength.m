function hours = daylength(declin, latitude)
% omega = daylight(declin, latitude)
%
% length of day, does not consider refraction
%
% Input
%   declin - solar declination, + in northen hemisphere
%   latitude - in degrees, + in northern hemisphere
%   (if both arguments are not scalars, they need to be same size)
%
% Output
%   hours - length of daylight

if ~isscalar(declin) && ~isscalar(latitude)
    assert(isequal(size(declin),size(latitude)),...
        'if not scalars, arguments must be of same size')
end

arg = -tand(latitude).*tand(declin);
% reset argument for total darkness (>1) or total sunlight(<-1)
t = abs(arg)>1;
if nnz(t)
    arg(t) = sign(arg(t));
end
hours = 2*acosd(arg)/15;

end