function [info] = singeotiffwrite(filename,X,print_cmap,R)
%SINGEOTIFF print a sinusodial projection geotiff
%   MODIS sinusoidal projection is not a predefined projected coordinate
%   system in the geotiff key system - this function manually lables all the
%   required keys correctly to print a MODIS tile to a geotiff. 
%  User-define projected coordinate systems may be defined by defining the
%Geographic Coordinate System, the coordinate transformation method and its
%associated parameters, as well as the planar system's linear units.
%
%
% Written by: Timbo Stillinger | 7/29/2020 | tcs@ucsb.edu

%Geographic Coordiante System info
% requires a user defined datum - MODIS sinusodial projection is calculated
% from a sphere that has the same area as the WGS 84 ellipsoid, which is
% a sphere of radius 6371007.181m
key.GeographicTypeGeoKey=32767;% user defined GCS
key.GeogPrimeMeridianGeoKey=8901;% PM_Greenwich
key.GeogLinearUnitsGeoKey=9001;%meters
key.GeogAngularUnitsGeoKey=9102;% degrees
key.GeogSemiMajorAxisGeoKey=6371007.181;
key.GeogSemiMinorAxisGeoKey=6371007.181;

%projection coordinate system info
key.GTModelTypeGeoKey = 1; % projected coordiante system
key.ProjectedCSTypeGeoKey = 32767; %user-defined Projected Coordinate System
key.GTRasterTypeGeoKey = 1; %pixel is area
key.ProjCoordTransGeoKey=24; % %siusodial projection is transofmation key 24
key.ProjFalseEastingGeoKey=0.0;
key.ProjFalseNorthingGeoKey=0.0;
key.ProjLinearUnitsGeoKey=9001; %meters
key.ProjCenterLongGeoKey=0;
key.ProjCenterLatGeoKey=0;


% geotiff can print data in pretty colors, if desired.
%could also scale data to save space here
if isempty(print_cmap)
    geotiffwrite(filename,X,R,'GeoKeyDirectoryTag',key)
else
    %scale to uint8
    scaledX=uint8(255*X);
    geotiffwrite(filename,scaledX,print_cmap,R,'GeoKeyDirectoryTag',key)
end

%check
info=geotiffinfo(filename);
end