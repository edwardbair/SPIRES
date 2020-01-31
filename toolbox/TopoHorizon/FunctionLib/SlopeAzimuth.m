function [ S ] = SlopeAzimuth(ZS)
% [ SlopeStruct ] = SlopeAzimuth(Zstruct)
%
%computes slopes and aspects on an elevation grid
%output in degrees, slopes upward from horizontal, aspects from south (+ccw)
%
%INPUT -
%   Zstruct - elevation structure (see GridStructure)
%
%OUTPUT
%   S - structure of slopes and azimuths, unsigned or signed integers in hundredths
%           of degrees, slopes upward from horizontal, aspects from south (+ccw)
%   S also includes some of the input information from the elevation
%       structure -- gridtype, RefMatrix, etc, along with
%       scale - factor to multipy slopes & aspects to convert to floating point numbers

if isfield(ZS,'scale')
    Z = double(ZS.Z)*ZS.scale;
else
    Z = ZS.Z;
end

switch ZS.gridtype
    case 'projected'
        [ slope,aspect ] = SlopeAzmProjected(ZS.RefMatrix, ZS.Projection, Z);
    case 'geographic'
        [slope, aspect ] = SlopeAzmGeographic(ZS.RefMatrix, Z, ZS.Geoid);
    case 'geolocated'
        [ slope,aspect ] = SlopeAzmGeolocated( ZS.Lat, ZS.Lon, Z, ZS.Geoid);
    otherwise
        error('elevation gridtype must be ''projected'', ''geographic'', or ''geolocated''')
end

% convert the slopes and aspects to single precision
S = ZS;
S = rmfield(S,'Z');
S.Slope = single(slope);
S.Aspect = single(aspect);

end