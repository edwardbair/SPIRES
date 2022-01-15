function proj = CaliforniaAlbers()
% proj = CaliforniaAlbers()
% if MATLAB version R2020b or later, return the projcrs object for the
% California Albers projection, also called Teale
% if MATLAB version R2020a or earlier, return the projection structure for
% the projection
%   This is an equaconic projection centered on the state of California and
%   is used by the California Data Exchange Center

if verLessThan('map','5.0')
    x = defaultm('eqaconicstd');
    x.mapparallels = [34 40.5];
    x.origin = [0 -120 0];
    x.falsenorthing = -4000000;
    E = referenceEllipsoid(7019);
    x.geoid = [E.SemimajorAxis E.Eccentricity];
    proj = defaultm(x);
else
    ESPG_code = 6414;
    proj = projcrs(ESPG_code);
end

end