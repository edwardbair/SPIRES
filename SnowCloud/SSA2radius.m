function radius = SSA2radius(ssa,radUnit)
% radius = SSA2radius(ssa,radUnit)
%converts the surface specific area (SSA) in m^2/kg to equivalent spherical
%radius
%
%Input
% ssa - specific surface area of snow grain, in m^2/kg
% radUnit - unit for output radius, typically 'mum' ('um') or 'mm' but any
% metric length unit is okay, see convertLengthUnits.m for all possibilities
%
%Output
% radius - equivalent spherical radius, in the units specified

rho = 917; % kg/m^3, for ice
radius = 3./(rho*ssa); % radius in meters
if ~strcmpi(radUnit,'m')
    radius = convertLengthUnits(radius,'m',radUnit);
end

end