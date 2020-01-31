function [ enthalpy ] = EnthalpyVaporization( T, substance )
% [ enthalpy ] = EnthalpyVaporization( T, 'water' or 'ice' )
%
% input
%   T - deg C, can be a vector or matrix
%   substance - either 'ice' or 'water'
%
% output
%   enthalpy, of vaporization in J/kg
%   (set to NaN where substance is 'ice' and T>0;
%
% Polynomial curve fits to Table 2.1. R. R. Rogers & M. K. Yau (1989),
% A Short Course in Cloud Physics 


switch substance
    case 'water'
        enthalpy = 1.e3*(2500.8-2.36*T + 0.0016*T.^2 - 0.00006*T.^3);
    case 'ice'
        enthalpy = 1.e3*(2834.1 - 0.29*T - 0.004*T.^2);
        neg = T>0;
        enthalpy(neg) = NaN;
    otherwise
        error('substance (2nd argument) must be either ''ice'' or ''water''')
end

end