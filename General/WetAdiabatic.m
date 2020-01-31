function [ salr ] = WetAdiabatic( T )
%WetAdiabatic saturated adiabatic lapse rate
% UNDER DEVELOPMENT!!
%
% Input
%   T - Kelvin, temperature of saturated air
%
% Output
%   salr - saturated adiabatic lapse rate, negative, in deg/km

g = 9.8076; % gravity
Hv = 2.5e6; % heat of vaporization
R = 8.314; % J/kmol/K
Mw = 18.015;
Md = 28.964;
Rd = R/Md;
Rw = R/Mw;
reps = 0.622;
Cp = 1005;
vp = SaturationVaporPressure(T,'water');
end