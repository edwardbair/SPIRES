function [ T ] = DewPoint( vp, substance, units )
% [ T ] = DewPoint( vp, substance, units )
%Dew point as a function of vapor pressure over water or ice
%
% input
%   vp - vapor pressure, can be a vector or matrix
%   substance - either 'ice' or 'water'
%   units - 'Pa', 'kPa', 'hPa' or 'mb' or 'mbar' (these three are the
%   same), 'bar'
%
% output
%   T - dew point temperature, Kelvin
%   (set to NaN where substance is 'ice' and vp > 611.658 Pa)

assert(strcmpi(substance,'ice') || strcmpi(substance,'water'),...
    'substance (2nd argument) must be either ''ice'' or ''water''')
assert(~any(vp(:)<0), 'vapor pressures must be non-negative')

T0 = 273.16;
es0 = 611.658; % sat vp in Pa at 273.16K

switch units
    case {'Pa','pa'}
        conversion = 1;
    case {'mb', 'mbar', 'hPa', 'hpa'}
        conversion = 100;
    case {'kPa','kpa'}
        conversion = 1000;
    case 'bar'
        conversion = 1.e5;
    otherwise
        error('units ''%s'' not recognized',units)
end

% convert to Pa
vp = conversion*vp;

if strcmpi(substance,'ice')
    T = nan(size(vp));
    z = vp==0;
    T(z) = 0;
    for k=1:numel(vp)
        if vp(k)<=es0 && vp(k)>0
            vpx = vp(k);
            try
                T(k) = fzero(@F_ice,[0 T0]);
            catch
                T(k) = fzero(@F_water,T0);
            end
        end
    end
else
    T = zeros(size(vp));
    for k=1:length(vp)
        if vp(k)>0
            vpx = vp(k);
            T(k) = fzero(@F_water,T0);
        end
    end
end

if strcmpi(substance,'ice') && any(vp(:)>es0)
    notIce = vp>es0;
    warning([num2str(nnz(notIce)) ' vapor pressure(s) > '...
        num2str(es0) ' Pa, calculating for water not ice'])
    T(notIce) = DewPoint(vp(notIce),'water','Pa');
end

    function vpw = F_water(T)
        vpw = SaturationVaporPressure(T,'water','Pa')-vpx;
    end

    function vpi = F_ice(T)
        vpi = SaturationVaporPressure(T,'ice','Pa')-vpx;
    end

end