function T = tauSnow(r,swe,Qext)
% T = tauSnow(r,swe,Qext)
% r in m
% swe in kg/m^2 (same as mm of swe)
% Qext Mie extinction efficiency
% (inputs that are not scalars must be same size)

rhoIce = 917; % kg/m^2
[r,swe,Qext] = checkSizes(r,swe,Qext);
T = (3*swe.*Qext)./(4*r*rhoIce);

end