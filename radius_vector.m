function rv = radius_vector( doy )
%radius_vector: approximation for Earth-Sun distance in AU
%
% Input
%   doy - calendar day-of year (can be a vector of days)
%
% Output
%   radius vector

a0 = 1.000128376223074;
a1 = -0.016662032168947;
b1 = -0.001013122469254;
a2 = -1.2783185533524e-4;
b2 = -2.152377408901758e-5;
w = 0.017191017439612;

rv = a0 + a1*cos(doy*w) + b1*sin(doy*w) + a2*cos(2*doy*w) + b2*sin(2*doy*w);

end

