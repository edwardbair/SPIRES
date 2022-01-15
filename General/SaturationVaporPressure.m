function [ vp, varargout ] = SaturationVaporPressure(T,substance,units)
% [ svp [,dsvp] ] = SaturationVaporPressure( T, substance, units )
%saturation vapor pressure over water or ice
%optionally, also derivative of saturation vapor pressure
%
% input
%   T - Kelvin temperature, can be a vector or matrix
%   substance - either 'ice' or 'water'
%   units - 'Pa', 'kPa', 'hPa' or 'mb' or 'mbar' (these three are the
%   same), 'bar'
%
% output
%   svp - saturation vapor pressure(s) in appropriate units corresponding to
%       temperature(s) (set to NaN where substance is 'ice' and T>273.15
% optional output
%   dsvp - 1st derivative of saturation vapor pressure(s) to temperature
%
%
% Equations from Hardy, B. (1998). ITS-90 formulations for vapor pressure,
% frostpoint temperature, dewpoint temperature, and enhancement factors in
% the range -100 TO +100 C. Paper presented at the Proceedings of the Third
% International Symposium on Humidity & Moisture, London.
% http://www.rhs.com/papers/its90form.pdf

assert(strcmpi(substance,'ice') || strcmpi(substance,'water'),...
    'substance (2nd argument) must be either ''ice'' or ''water''')
assert(~all(isnan(T(:))),'all temperatures are NaN')
assert(~any(T(:)<0 & ~isnan(T(:))), 'Temperatures must be in Kelvin')

persistent T0 k0 k1 k2 k3 k4 k5 g0 g1 g2 g3 g4 g5 g6 g7

if isempty(T0)
    T0 = 273.16;
    k0 = -5.8666426e3;
    k1 = 2.232870244e1;
    k2 = 1.39387003e-2;
    k3 = -3.4262402e-5;
    k4 = 2.7040955e-8;
    k5 = 6.7063522e-1;
    g0 = -2.8365744e3;
    g1 = -6.028076559e3;
    g2 = 1.954263612e1;
    g3 = -2.737830188e-2;
    g4 = 1.6261698e-5;
    g5 = 7.0229056e-10;
    g6 = -1.8680009e-13;
    g7 = 2.7150305;
end

switch substance
    case 'ice'
        loges = (k0+k1*T+T.^2.*(k2+T.*(k3+k4.*T)))./T + k5*log(T);
    case 'water'
        loges = (g0+g1*T+T.^2.*(g2+T.*(g3+T.*(g4+T.*(g5+g6.*T)))))./T.^2 +...
            g7*log(T);
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

vp = conversion * exp(loges);
vp(T==0) = 0;

if nargout>1
    assert(nargout==2,'only one optional output argument')
    switch substance
        case 'ice'
            dsvp = conversion*exp(loges).*...
                (-k0+k5*T+T.^2.*(k2+T.*(2*k3+3*k4.*T)))./T.^2;
        case 'water'
            dsvp = conversion*exp(loges).*...
                (g7./T+(-2*g0-g1*T+T.^3.*(g3+T.*(2*g4+T.*(3*g5+4*g6.*T))))./T.^3);
    end
    varargout{1} = dsvp;
end

if strcmpi(substance,'ice') && any(T(:)>T0)
    notIce = T>T0;
    warning([num2str(nnz(notIce)) ' temperature(s) > ' num2str(T0)...
        ' treated as water not ice'])
    if nargout>1
        [vp(notIce),dsvp(notIce)] = SaturationVaporPressure(T(notIce),'water',units);
        varargout{1} = dsvp;
    else
        vp(notIce) = SaturationVaporPressure(T(notIce),'water',units);
    end
end

end