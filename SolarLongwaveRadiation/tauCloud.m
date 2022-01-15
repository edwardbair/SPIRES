function T = tauCloud(r,waterEquiv,Qext)
% T = tauCloud(r,waterEquiv,Qext)
% r in m
% waterEquiv in kg/m^2, same as mm

rhoWater = 1000; % kg/m^3
[r,waterEquiv,Qext] = checkSizes(r,waterEquiv,Qext);
T = (3*waterEquiv.*Qext)./(4*r*rhoWater);

end