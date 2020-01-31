function [ vp ] = SaturationVaporPressureBA( T, substance, units )
% [ vp ] = SaturationVaporPressure( T, substance, units )
%saturation vapor pressure over water or ice
%
% input
%   T - Kelvin temperature, can be a vector or matrix
%   substance - either 'ice' or 'water'
%   units - 'Pa', 'kPa', 'hPa' or 'mb' or 'mbar' (these three are the
%   same), 'bar'
%
% output
%   vp - saturation vapor pressure(s) in appropriate units corresponding to
%       temperature(s) (set to NaN where substance is 'ice' and T>273.15
%
% Old version equations are from Bohren, C. F., and B. A. Albrecht (1998),
% Atmospheric Thermodynamics, Oxford University Press.
%
% New version from Hardy, B. (1998). ITS-90 formulations for vapor pressure,
% frostpoint temperature, dewpoint temperature, and enhancement factors in
% the range -100 TO +100 C. Paper presented at the Proceedings of the Third
% International Symposium on Humidity & Moisture, London.
% http://www.rhs.com/papers/its90form.pdf



assert(strcmpi(substance,'ice') || strcmpi(substance,'water'),...
    'substance (2nd argument) must be either ''ice'' or ''water''')
assert(nnz(T(T<0)) == 0, 'Temperatures must be in Kelvin')

T0 = 273.15;
es0 = 611; % sat vp in Pa at 273.15K
switch substance
    case 'ice'
        logs = 6293.*(1/T0 - 1./T) - 0.555*log(T/T0);
        ok = T <= T0;
        if nnz(~ok)>0
            logs(~ok) = NaN;
        end
    case 'water'
        logs = 6808.*(1/T0 - 1./T) - 5.09*log(T/T0);
    otherwise
        error('substance (2nd argument) must be either ''ice'' or ''water''')
end

switch units
    case {'Pa','pa'}
        conversion = 1;
    case {'mb', 'mbar', 'hPa', 'hpa'}
        conversion = 1/100;
    case {'kPa','kpa'}
        conversion = 1/1000;
    case 'bar'
        conversion = 1.e-5;
    otherwise
        error('units ''%s'' not recognized',units)
end

vp = conversion * es0 * exp(logs);

end