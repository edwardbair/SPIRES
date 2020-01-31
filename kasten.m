function airmass = kasten( mu0 )
% airmass = kasten( mu0 )
%Kasten Relative optical airmass from cosine of the solar zenith angle
%   Equation from Kasten, F., and A. T. Young (1989), Revised optical air
%   mass tables and approximation formula, Applied Optics, 28, 4735-4738,
%   doi: 10.1364/AO.28.004735.

% coefficients
a = 0.50572;
b = 6.07995;
c = 1.6364;

% solar elevation in degrees, from horizon upward
gam = asind(mu0);

% Kasten-Young equation - set to NaN if below horizon
t = mu0 < 0;
if nnz(t)
    airmass = NaN(size(mu0));
    airmass(~t) = 1 ./(sind(gam(~t))+a.*(gam(~t)+b).^(-c));
else
    airmass = 1 ./(sind(gam)+a.*(gam+b).^(-c));
end

% slight correction for overhead sun
airmass(airmass<1) = 1;

end