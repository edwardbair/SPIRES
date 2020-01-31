function [S ,I] = read_AVNG_ENVIhdr(hdrFile)
%read the hdr file for the AVIRIS-NG data in ENVI format and return
%information
%Output
%   S - header file for analysis of AVIRIS-NG files
%   I - output from envihdrread

%descriptions in the headers themselves are arcane, so use these based on
%the filename
highLevelDescription = {'corr','h2o','rdn'}; % corrected values, or radiances
lowLevelDescription = {'img','obs','igm','loc'};
wordDescription = {'AVIRIS-NG atmospherically corrected reflectance',...
    'AVIRIS-NG estimated H2O vapor, liquid, ice in cm',...
    'AVIRIS-NG measured radiance in uW nm-1 cm-2 sr-1',...
    'AVIRIS-NG lengths and angles',...
    'AVIRIS-NG projection coordinates',...
    'AVIRIS-NG geographic coordinates'}';
for n=1:length(highLevelDescription)
    if contains(hdrFile,highLevelDescription{n})
        switch n
            case 1
                S.description = wordDescription{1};
                S.outputVariableName = 'reflectance';
                S.badNegative = true;
                S.outerLeave = 'BIP';
                break
            case 2
                S.description = wordDescription{2};
                S.outputVariableName = 'h2o';
                S.badNegative = true;
                S.outerLeave = 'BSQ';
                break
            case 3
                for m=1:length(lowLevelDescription)
                    if contains(hdrFile,lowLevelDescription{m})
                        switch m
                            case 1
                                S.description = wordDescription{3};
                                S.outputVariableName = 'radiance';
                                S.badNegative = true;
                                S.outerLeave = 'BIP';
                            case 2
                                S.description = wordDescription{4};
                                S.outputVariableName = 'obs';
                                S.badNegative = false;
                                S.outerLeave = 'BSQ';
                            case 3
                                S.description = wordDescription{5};
                                S.outputVariableName = 'projCoord';
                                S.badNegative = false;
                                S.outerLeave = 'BSQ';
                            case 4
                                S.description = wordDescription{6};
                                S.outputVariableName = 'geogCoord';
                                S.badNegative = false;
                                S.outerLeave = 'BSQ';
                        end
                    end
                end
                break;
        end
    end
end

I = envihdrread(hdrFile);

% samples, lines, bands
S.samples = I.samples;
S.lines = I.lines;
S.bands = I.bands;
S.interleave = I.interleave;

% data type
iscx = false;
switch I.data_type
    case 1
        S.data_type = 'uint8';
    case 2
        S.data_type= 'int16';
        case3
        S.data_type= 'int32';
    case 4
        S.data_type= 'single';
    case 5
        S.data_type= 'double';
    case 6
        iscx=true;
        S.data_type= 'single';
    case 9
        iscx=true;
        S.data_type= 'double';
    case 12
        S.data_type= 'uint16';
    case 13
        S.data_type= 'uint32';
    case 14
        S.data_type= 'int64';
    case 15
        S.data_type= 'uint64';
    otherwise
        error(['File type number: ',num2str(I.data_type),' not supported']);
end
assert(~iscx,'complex data type not supported (yet)')

% bad values
if contains(S.data_type,'single','IgnoreCase',true) ||...
        contains(S.data_type,'double','IgnoreCase',true)
    S.ignoreValue = NaN;
elseif contains(S.data_type,'uint','IgnoreCase',true)
    S.ignoreValue = intmax(S.data_type);
else
    S.ignoreValue = intmin(S.data_type);
end

% wavelength
if isfield(I,'wavelength')
    line = I.wavelength;
    line(line=='{' | line=='}') = [];
    line=textscan(line,'%s','Delimiter',',');
    line=line{:};
    line=strtrim(line);
    S.wavelength = str2double(line);
end

% referencing matrix
if isfield(I,'map_info')
    M = I.map_info;
    M.cornerX = I.map_info.mapx;
    M.cornerY = I.map_info.mapy;
    M.pixelSize = [I.map_info.dx I.map_info.dy];
    if isfield(I.map_info,'rotation')
        M.rotation = I.map_info.rotation;
    else
        M.rotation = 0;
    end
    
    % projection information -- WARNING: specific to the AVIRIS-NG Indian data
    if isfield(M,'projection')
        if contains(M.projection,'Albers')
            S.projection = HimachalPradeshAlbers;
        elseif contains(M.projection,'utm','IgnoreCase',true)
            S.projection = HimachalPradeshUTM;
        else
            warning('projection not recognized')
            S.projection = [];
        end
    else
        S.projection = [];
    end
    mfield = {'cornerX','cornerY','pixelSize','rotation'};
    for k=1:length(mfield)
        S.mapping.(mfield{k}) = M.(mfield{k});
    end
    S.doMap = true;
else
    S.doMap = false;
end

if contains(hdrFile,'h2o')
    S.bandName = {'Column water vapor','Liquid absorption path','Ice absorption path'};
elseif contains(hdrFile,'obs')
    S.bandName = {'Path length (m)','To-sensor azimuth (0 to 360 degrees cw from N)',...
        'To-sensor zenith (0 to 90 degrees from zenith)',...
        'To-sun azimuth (0 to 360 degrees cw from N)',...
        'To-sun zenith (0 to 90 degrees from zenith)','Solar phase',...
        'Slope','Aspect','Cosine(i)','UTC Time Offset','Earth-sun distance (AU)'};
    S.badNegative = true(1,length(S.bandName));
    S.badNegative(end-3:end-2) = false;
elseif contains(hdrFile,'igm')
    S.bandName = {'Easting','Northing','Elevation'};
elseif contains(hdrFile,'loc')
    S.bandName = {'Longitude','Latitude','Elevation'};
end

% bad value
if isfield(I,'data_ignore_value')
    S.badValue = I.data_ignore_value;
end

end