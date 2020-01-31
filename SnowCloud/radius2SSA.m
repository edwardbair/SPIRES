function ssa = radius2SSA(radius,radUnit)
% ssa = radius2SSA(radius,radUnit)
%converts equivalent spherical radius to SSA (specific surface area) in
%m^2/kg
%
%Input
% radius - equivalent spherical radius, in the units specified
% radUnit - unit for output radius, typically 'mum' ('um') or 'mm' but any
%   metric length unit is okay, see convertUnits.m for all possibilities
%
%Output
% ssa - specific surface area of snow grain, in m^2/kg


rho = 917; % kg/m^3, for ice
if ~strcmpi(radUnit,'m')
    radius = convertLengthUnits(radius,radUnit); % to meters, default
end
ssa = 3./(rho*radius); % radius in meters

end