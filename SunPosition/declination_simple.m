% approximate (better than 1/10th degree) solar declination in degrees
function dec=declination_simple(doy)
% dec=declination_simple(doy)
% input doy numerical day of year, starting 1 at 01 January
A0=0.292033450096867;
A=[-22.982158535385800 -0.287744385009193 -0.134810057611257]';
B=[3.751757849830815 0.063793006836552 0.103899662868187]';
omega=0.017115168775574;

angles=doy*omega;
dec=A0+A(1)*cos(angles)+B(1)*sin(angles)+A(2)*cos(2*angles)+...
    B(2)*sin(2*angles)+A(3)*cos(3*angles)+B(3)*sin(3*angles);
end