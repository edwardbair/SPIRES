function proj = MODISsinusoidal()
% proj = MODISsinusoidal()
% if MATLAB version R2020b or later, return the projcrs object for the
% sinusoidal projection used by MODIS and other NASA satellite sensors
% (uses the sphere of the "authalic" radius corresponding to WGS 1984)
% if MATLAB version R2020a or earlier, return the projection structure for
% the projection

if verLessThan('map','5.0')
    MODISgeoid = [6.371007181e+06 0];
    mstruct = defaultm('sinusoid');
    mstruct.geoid = MODISgeoid;
    proj = defaultm(mstruct);
else
    % (MATLAB could use a more explicit method to create a custom projcrs)
    % (the following gives the same result as the mstruct version)
    wkt = "PROJCS[""MODIS Sinusoidal"",BASEGEOGCRS[""User"",DATUM[""World Geodetic Survey 1984"",SPHEROID[""Authalic_Spheroid"",6371007.181,0.0]],PRIMEM[""Greenwich"",0.0],UNIT[""Degree"",0.0174532925199433]],PROJECTION[""Sinusoidal""],PARAMETER[""False_Easting"",0.0],PARAMETER[""False_Northing"",0.0],PARAMETER[""Central_Meridian"",0.0],UNIT[""Meter"",1.0]]";
    % create the projcrs from the edited wkt
    proj = projcrs(wkt);
end
end