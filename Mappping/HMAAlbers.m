function proj = HMAAlbers()
% proj = HMAAlbers()
% if MATLAB version R2020b or later, return the projcrs object for David
% Shean's custom Albers equal area projection for High Mountain Asia
% if MATLAB version R2020a or earlier, return the projection structure for
% the projection
% Documentation at https://nsidc.org/data/highmountainasia

if verLessThan('map','5.0')
    p = defaultm('eqaconicstd');
    p.mapparallels = [25 47];
    p.origin = [36 85 0];
    E = referenceEllipsoid(7030);
    p.falsenorthing = 0;
    p.falseeasting = 0;
    p.geoid = [E.SemimajorAxis E.Eccentricity];
    proj = defaultm(p);
else
    % hack hack: I started with the US Albers projection, got its wkt (well
    % known string), edited the wkt, and then use it to create the projcrs
    % (MATLAB could use a more explicit method to create a custom projcrs)
    wkt = "PROJCRS[""High Mountain Asia Albers Equal Area Conic"",BASEGEOGCRS[""WGS84"",DATUM[""World Geodetic Survey 1984"",ELLIPSOID[""WGS 1984"",6378137,298.257223563,LENGTHUNIT[""metre"",1]]],PRIMEM[""Greenwich"",0,ANGLEUNIT[""Degree"",0.0174532925199433]]],CONVERSION[""High Mountain Asia Albers Equal Area Conic"",METHOD[""Albers Equal Area"",ID[""EPSG"",9822]],PARAMETER[""Latitude of false origin"",36,ANGLEUNIT[""Degree"",0.0174532925199433],ID[""EPSG"",8821]],PARAMETER[""Longitude of false origin"",85,ANGLEUNIT[""Degree"",0.0174532925199433],ID[""EPSG"",8822]],PARAMETER[""Latitude of 1st standard parallel"",25,ANGLEUNIT[""Degree"",0.0174532925199433],ID[""EPSG"",8823]],PARAMETER[""Latitude of 2nd standard parallel"",47,ANGLEUNIT[""Degree"",0.0174532925199433],ID[""EPSG"",8824]],PARAMETER[""Easting at false origin"",0,LENGTHUNIT[""metre"",1],ID[""EPSG"",8826]],PARAMETER[""Northing at false origin"",0,LENGTHUNIT[""metre"",1],ID[""EPSG"",8827]]],CS[Cartesian,2],AXIS[""(E)"",east,ORDER[1],LENGTHUNIT[""metre"",1]],AXIS[""(N)"",north,ORDER[2],LENGTHUNIT[""metre"",1]],USAGE[SCOPE[""unknown""],AREA[""High Mountain Asia""],BBOX[23,60,54,110]]]";
    % create the projcrs from the edited wkt
    proj = projcrs(wkt);
end
end