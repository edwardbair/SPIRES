function [ RefMatrix, varargout ] = sinusoidProjMODtile( tile )
% [ RefMatrix ] = sinusoidProjMODtile( tile )
% [ RefMatrix,ProjectionStructure ] = sinusoidProjMODtile( tile )
% [ RefMatrix,ProjectionStructure,RasterReference] = sinusoidProjMODtile( tile )
% [ RefMatrix,ProjectionStructure,RasterReference,GeoKeyDirectoryTag] = sinusoidProjMODtile( tile )
% [ RefMatrix,ProjectionStructure,RasterReference,GeoKeyDirectoryTag, latlonCorners ] = sinusoidProjMODtile( tile )
%
%Coordinate and projection information for any MODIS tile
%
% Input
%   tile - 6-character MODIS tile designation in form 'hNNvNN'
%
% Output
%   RefMatrix - structure containing referencing matrices for 250 m, 500 m,
%       and 1 km MODIS data
% Optional output, if specified in output arguments as above
%   ProjectionStructure - sinusoidal MATLAB projection structure to use in
%       converting from projection coordinates to latitude and longitude
%   RasterReference - structure containing the mapRasterReference objects
%       for each of the MODIS resolutions
%   crs now appended as of 2021a to RR
%   GeoKeyDirectoryTag - structure containing the various GeoKeys for input
%       to a GeoTIFF file

%   latlonCorners - latitudes and longitudes for each corner, counter-
%       clockwise from upper left, 4x2 matrix with latitudes in column 1,
%       longitudes in column 2
%
% The function assumes that the grid origin is in the upper left

persistent already mstruct

% tile sizes for 1km, 500m, 250m, 375m(VIIRS) - hardwired
tilesize = [1200 1200; 2400 2400; 4800 4800; 3000 3000];

% loop through the possible arguments
if nargout>1
    varargout = cell(nargout-1,1);
end
if iscell(tile)
    tile = char(tile);
end
for k=1:nargout
    switch k
        case 1 % referencing matrix, required
            % x- and y-coordinates of this tile's upper left corner
            [ULx,ULy,tileHeight,tileWidth] = MODtile2xy(tile);
            pixelHeight = tileHeight./tilesize(:,1);
            pixelWidth = tileWidth./tilesize(:,2);
            
            % coordinates of center of upper left pixel (these are vectors
            % because pixelWidth and pixelHeight are vectors)
            x11 = ULx+pixelWidth/2;
            y11 = ULy-pixelHeight/2; % negative because y(n)>y(n+1)
            
            % Referencing matrices
            refmat = {'RefMatrix_1km','RefMatrix_500m','RefMatrix_250m', 'RefMatrix_375m'};
            assert(length(refmat)==length(pixelWidth),...
                'number of referencing matrices (%d) does not equal number of tile sizes (%d)',...
                length(refmat),length(pixelWidth))
            for gn=1:length(refmat)
                RefMatrix.(refmat{gn}) = makerefmat(x11(gn),y11(gn),...
                    pixelWidth(gn),-pixelHeight(gn));
            end
        case 2 % projection structure
            % load the projection structure first time we need it
            if isempty(already)
                already = true;
                % hardwired, authalic radius corresponding to WGS84
                MODISgeoid = [6.371007181e+06 0];
                mstruct = defaultm('sinusoid');
                mstruct.geoid = MODISgeoid;
                mstruct = defaultm(mstruct);
            end
            varargout{k-1} = mstruct;
        case 3 % map raster references
            rasterref = {'RasterReference_1km','RasterReference_500m','RasterReference_250m', 'RasterReference_375m'};
            assert(length(rasterref)==length(refmat),...
                'bug in code, length(rasterref)=%d, length(refmat)=%d',...
                length(rasterref),length(refmat))
            %R2021 
            crs=projcrs(['PROJCS["MODIS Sinusoidal",BASEGEOGCRS["User",',...
                'DATUM["World Geodetic Survey 1984",',...
                'SPHEROID["Authalic_Spheroid",6371007.181,0.0]],',...
                'PRIMEM["Greenwich",0.0],',...
                'UNIT["Degree",0.0174532925199433]],',...
                'PROJECTION["Sinusoidal"],',...
                'PARAMETER["False_Easting",0.0],',...
                'PARAMETER["False_Northing",0.0],',...
                'PARAMETER["Central_Meridian",0.0],UNIT["Meter",1.0]]']);

            for gn=1:length(refmat)
                RasterReference.(rasterref{gn}) =...
                    refmatToMapRasterReference(RefMatrix.(refmat{gn}),...
                    tilesize(gn,:));
                RasterReference.(rasterref{gn}).ProjectedCRS=crs;
           
            end
            
            
            varargout{k-1} = RasterReference;
        case 4 % GeoKeyDirectoryTag
            CT_Sinusoidal = 24;
            ModelTypeProjected = 1;
            RasterPixelIsArea = 1;
            LinearMeter = 9001;
            AngularDegree = 9102;
            key.GTModelTypeGeoKey = ModelTypeProjected;
            key.GTRasterTypeGeoKey = RasterPixelIsArea;
            key.GeogLinearUnitsGeoKey = LinearMeter;
            key.GeogAngularUnitsGeoKey = AngularDegree;
            key.ProjectedCSTypeGeoKey = CT_Sinusoidal;
            varargout{k-1} = key;
        case 5 % lat-lon of corners
            [clat,clon] = MODtile2latlon(tile);
            varargout{k-1} = [clat clon];

        otherwise
            error('too many output arguments')
    end
end
end