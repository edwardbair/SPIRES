function S = anglesMOD1B(file)
% S = anglesMOD1B(file)
%sensor and solar angles, and latitude and longitude from MODIS Level 1B file
%
%Input
%   file - MODIS filename
%
%Output (angles in degrees)
%   S - structure with elements
%       CosineSensorZenith,
%       SensorAzimuth,
%       CosineSolarZenith,
%       SolarAzimuth,
%       Latitude,
%       Longitude
% Output azimuths are from SOUTH, positive ccw

loc = {'Latitude','Longitude'};
ang = {'SensorZenith','SensorAzimuth','SolarZenith','SolarAzimuth'};
angParse = {'CosineSensorZ','SensorAzimuthS','CosineSolarZ','SolarAzimuthS'};

% find valid ranges, scales and offsets
I = hdfinfo(file,'hdf');
V = I.Vgroup.Vgroup;
for k=1:length(V)
    SDS = V(k).SDS;
    for n=1:length(SDS)
        A = SDS(n).Attributes;
        switch SDS(n).Name
            case loc
                for m=1:length(A)
                    if strcmp(A(m).Name,'valid_range')
                        Range.(SDS(n).Name) = double(A(m).Value);
                    end
                end
            case ang
                for m=1:length(A)
                    switch A(m).Name
                        case 'scale_factor'
                            Scale.(SDS(n).Name) = double(A(m).Value);
                        case 'valid_range'
                            Range.(SDS(n).Name) = double(A(m).Value);
                    end
                end
        end
    end
end

% read values and set NaN where out of range
for k=1:length(loc) % lat/lon are floats
    assert(isfield(Range,loc{k}),'range for %s not found',loc{k})
    X = hdfread(file,loc{k});
    X(X<min(Range.(loc{k}))) = NaN;
    X(X>max(Range.(loc{k}))) = NaN;
    S.(loc{k}) = X;
end
for k=1:length(ang) % angles are scaled integers
    X = hdfread(file,ang{k});
    thisAng = single(X)*Scale.(ang{k});
    thisAng(X<min(Range.(ang{k}))) = NaN;
    thisAng(X>max(Range.(ang{k}))) = NaN;
    if strfind(ang{k},'Zenith') % cosine
        cX = cosd(thisAng);
        cX(cX<0) = 0;
        S.(angParse{k}) = cX;
    elseif strfind(ang{k},'Azimuth') % change to ccw from S
        S.(angParse{k}) = 180-thisAng;
        S.(angParse{k})(cX<0) = NaN;
    else
        error('ang %s not recognized',ang{k})
    end
end
end