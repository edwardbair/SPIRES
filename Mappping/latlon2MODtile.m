function [ MODIStile ] = latlon2MODtile( latitude,longitude )
% [ MODIStile ] = latlon2MODtile( latitude,longitude )
%returns MODIS tiles, which generally overlap, from latitude and longitude
%
% Input
%   latitude and longitude in degrees, either scalars or vectors of equal
%   length
%       latitudes are negative in the southern hemisphere
%       longitudes are negative in the western hemisphere
%
% Output
%   MODIStile
%       as a character string in form 'hNNvNN' if latitude and
%       longtiude are scalars
%       as a cell vector if latitude and longtitude are vectors

persistent already mstruct
if isempty(already)
    already = true;
    % hardwired, authalic radius corresponding to WGS84
    MODISgeoid = [6.371007181e+06 0];
    mstruct = defaultm('sinusoid');
    mstruct.geoid = MODISgeoid;
    mstruct = defaultm(mstruct);
end

% preliminaries
EarthCircumference = 2*pi*earthRadius;
nTilesHorizontal = 36;
tileWidth = EarthCircumference/nTilesHorizontal;
tileHeight = tileWidth;
nTilesVertical = 18;

% x- and y-coordinates in sinusoidal projection
assert(isequal(size(latitude),size(longitude)),...
    'vectors of latitude and longitude must be same size')
% convert to column vectors
if length(latitude)>1 && isrow(latitude)
    latitude = latitude';
    longitude = longitude';
end
[x,y] = mfwdtran(mstruct,latitude,longitude);

% h and v coordinates, floating point (inverse of equations in MODtile2xy)
h = (EarthCircumference+2*x)/(2*tileWidth);
v = -(EarthCircumference-4*nTilesVertical*tileHeight+4*y)/(4*tileHeight);

if isscalar(h)
    MODIStile = ['h', num2str(floor(h),'%02d')...
        'v' num2str(floor(v),'%02d')];
else
    MODIStile = cell(length(h),1);
    for c=1:length(h)
        MODIStile{c} = ['h', num2str(floor(h(c)),'%02d')...
            'v' num2str(floor(v(c)),'%02d')];
    end
    
end

end