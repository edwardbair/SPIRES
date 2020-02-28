function [IVT] = calculateIVT(q,u,v,g)
%calculate the IVT from 1000 hPa to 300 hPa
% global data: ERA5 hourly data on pressure levels from 1979 to present
% CA data: West WRF CW3E
% q is the specific humidity in kg / kg
% u is the zonal wind in ms?1
% v is the meridional wind in ms?1
%g is the acceleration due to gravity.


%sum specific humidty*u for each pressure layer 300-1000
%pressure layers are unit spacing. Y = q*u vector for each pressure layer
qu=q.*u;
qv =q.*v;
intU = trapz(qu);
intV = trapz(qv);

IVT = sqrt((intU./g)^2+(intV./g)^2);
end