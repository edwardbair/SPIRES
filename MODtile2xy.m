function [ ULx,ULy,tileHeight,tileWidth ] = MODtile2xy( tile )
% [ ULx,ULy,tileHeight,tileWidth ] = MODtile2xy( tile )
%upper left coordinates of MODIS tile, and its height and width
%
% Input
%   tile - MODIS tile in form 'hNNvNN'
%
% Output
%   ULx, ULy - coordinates of upper left corner of tile in sinusoidal
%       projection (meters)
%   tileHeight, tileWidth - height and width of tile, meters
%
% Acknowledgment: Duke's Marine Spatial Geoecology Tools
% http://code.env.duke.edu/projects/mget/wiki/SinusoidalMODIS

% on first pass, make table of off-Earth tiles
persistent already T
if isempty(already)
    already = true;
    nbad = [28 22 18 12 8 4 2 0 0 0 0 2 4 8 12 18 22 28];
    firsth = 0;
    lasth = 35;
    v = [];
    h = [];
    for k=1:length(nbad)
        if nbad(k)>0
            badv = k-1;
            v=cat(1,v,zeros(nbad(k),1)+badv);
            badh = cat(1,(firsth:nbad(k)/2-1)',(1+lasth-nbad(k)/2:lasth)');
            h = cat(1,h,badh);
        end
    end
    T = table(v,h);
end

% authalic radius corresponding to WGS84 hardwired for MODIS sinusoidal
% projection
earthRadius = 6.371007181e+06;
nTilesHorizontal = 36;
nTilesVertical = 18;
EarthCircumference = 2*pi*earthRadius;
tileWidth = EarthCircumference/nTilesHorizontal;
tileHeight = tileWidth;

% tile coordinates
horzTile = str2double(tile(2:3));
vertTile = str2double(tile(5:6));
% check if off Earth
if ~isempty(find(vertTile==T.v & horzTile==T.h,1))
    error('tile %s is off Earth',tile)
end

% x,y of upper left corner
ULx = -EarthCircumference/2+horzTile*tileWidth;
ULy = -EarthCircumference/4+tileHeight*(nTilesVertical-vertTile);

end