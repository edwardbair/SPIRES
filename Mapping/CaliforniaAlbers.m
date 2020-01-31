function projstruct = CaliforniaAlbers()
% projstruct = CaliforniaAlbers()
%return the projection structure for the California Albers projection (also
%called Teale)
%   This is an equaconic projected centered on the state of California and
%   is used by the California Data Exchange Center

x = defaultm('eqaconicstd');
x.mapparallels = [34 40.5];
x.origin = [0 -120 0];
x.falsenorthing = -4000000;
E = wgs84Ellipsoid();
x.geoid = [E.SemimajorAxis E.Eccentricity];
projstruct = defaultm(x);

end