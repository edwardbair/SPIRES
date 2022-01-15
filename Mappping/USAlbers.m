function proj = USAlbers()
% proj = USAlbers()
% if MATLAB version R2020b or later, return the projcrs object for the
% US Albers projection, used for national maps that cover the whole CONUS
% if MATLAB version R2020a or earlier, return the projection structure for
% the projection

% note: The ESRI definition differs from the EPSG definition:
% ESRI (code 102003) puts the parallels at [29.5 45.5] and the origin at [37.5 -96]
%   and uses the NAD83 datum
% EPSG (code 6350) uses the same parallels but orign at [23 -96] and the
%   NAD83 datum
% USGS, for the ARD products, uses the EPSG definition but the WGS 84 datum

% mapping toolbox version 5.0, which has the projcrs function, was introduced in R2020b
if verLessThan('map','5.0')
    p = defaultm('eqaconicstd');
    p.mapparallels = [29.5 45.5];
    p.origin = [23 -96 0];
    p.falsenorthing = 0;
    p.falseeasting = 0;
    E = referenceEllipsoid(7019);
    p.geoid = [E.SemimajorAxis E.Eccentricity];
    proj = defaultm(p);
else
    wkt = "PROJCRS[""US Albers"",BASEGEOGCRS[""WGS 84"",DATUM[""World Geodetic System 1984"",ELLIPSOID[""WGS 84"",6378140,298.256999999996,LENGTHUNIT[""metre"",1]]],PRIMEM[""Greenwich"",0,ANGLEUNIT[""degree"",0.0174532925199433]],ID[""EPSG"",4326]],CONVERSION[""unnamed"",METHOD[""Albers Equal Area"",ID[""EPSG"",9822]],PARAMETER[""Latitude of false origin"",23,ANGLEUNIT[""degree"",0.0174532925199433],ID[""EPSG"",8821]],PARAMETER[""Longitude of false origin"",-96,ANGLEUNIT[""degree"",0.0174532925199433],ID[""EPSG"",8822]],PARAMETER[""Latitude of 1st standard parallel"",29.5,ANGLEUNIT[""degree"",0.0174532925199433],ID[""EPSG"",8823]],PARAMETER[""Latitude of 2nd standard parallel"",45.5,ANGLEUNIT[""degree"",0.0174532925199433],ID[""EPSG"",8824]],PARAMETER[""Easting at false origin"",0,LENGTHUNIT[""metre"",1],ID[""EPSG"",8826]],PARAMETER[""Northing at false origin"",0,LENGTHUNIT[""metre"",1],ID[""EPSG"",8827]]],CS[Cartesian,2],AXIS[""easting"",east,ORDER[1],LENGTHUNIT[""metre"",1,ID[""EPSG"",9001]]],AXIS[""northing"",north,ORDER[2],LENGTHUNIT[""metre"",1,ID[""EPSG"",9001]]]]";
    proj = projcrs(wkt);
end

end