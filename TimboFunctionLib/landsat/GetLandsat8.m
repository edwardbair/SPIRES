function [ X, outputFolder, varargout] = GetLandsat8( filepath, whichVariable, varargin )
%extract variable from Landsat8 input file - WORKS WITH LANDSAT 8
%COLLECTION 1 DATA
% [ X, outputFolder, varargout] = GetLandsat8_Collection1( filepath, whichVariable, varargin )
%
% Input
%   filepath - name of Landsat8 file or folder. works with .tif, .txt, .tar, or .tar.gz files
%   whichVariable - choose among (case insensitive)
%           'SolarZenith','SolarAzimuth', 'num_obs'
%           'band1', 'band2', 'band3', ... 'band11','allbands'
%           'BQA', 'ReferenceObject', 'metadata'
%
%
% Output
%   X - Variable scaled to appropriate values and units (the angles and obs
%       variables are returned as scalars,  BQA is
%       returned as structures, band values as raw digital numbers)
%   outputFolder - full path name to output folder that contains the
%       Landsat 8 files (created if filepath is a .tar or .tar.gz file,
%       otherwise same as filepath)
%
%
% optional input
%   masterOutputFolder - folder where new folder with name Lxxxxx will go
%       (option is folder where inputL8file lives)
%
%optional output
%   geotiffread return - can specify to return referenceing matrix w/
%   bounding box [REFMAT, BBOX] or referencing matrix [R]
%
%
% 8/16/14 - add ability to get thermal bands 10 & 11.
% 9/3/14 - add 'metadata' case to whichVariable - outputs all MTL file info
% 9/4/14 - add 'geotiffinfo' return case

numOut = 2;

if iscell(filepath)
    filepath = char(filepath);
end
assert(ischar(filepath),'input filepaths must be cell or character arrays');

%extract info from filepath
[pathstr, ~ , sceneID] = filepathInfo(filepath);

% put where it is already, unless specified
if nargin>2
    outputFolder = fullfile(varargin{1},sceneID);
else
    outputFolder = pathstr;
end

%find MTL filename and MTL metadata
[MTL_list, MTL_filename, dataFolder] = GetMTL(filepath, outputFolder);

%figure out which version of Landsat 8 data
if isfield(MTL_list.METADATA_FILE_INFO,'COLLECTION_NUMBER')
    datatype =['collection' num2str(MTL_list.METADATA_FILE_INFO.COLLECTION_NUMBER)];
else
    datatype = 'precollection';
end


%grab/compute requested variable from metadata or band file
switch lower(whichVariable)
    case lower('solarZenith')
        X = 90 - MTL_list.IMAGE_ATTRIBUTES.SUN_ELEVATION;
    case lower('solarAzimuth')
        X = MTL_list.IMAGE_ATTRIBUTES.SUN_AZIMUTH;
    case lower('metadata')
        X = MTL_list;
    case lower('GeoTIFFinfo_30m')
        %geotiffinfo strucutre for 30m grid bands gathered from band 1
        band = 'FILE_NAME_BAND_1';
        bandFile = fullfile(pathstr, MTL_list.PRODUCT_METADATA.(band));
        X = geotiffinfo(bandFile);
    case 'num_obs_30m'
        X = MTL_list.PRODUCT_METADATA.REFLECTIVE_LINES .* ...
            MTL_list.PRODUCT_METADATA.REFLECTIVE_SAMPLES;
    case 'num_obs_15m'
        X = MTL_list.PRODUCT_METADATA.PANCHROMATIC_LINES .* ...
            MTL_list.PRODUCT_METADATA.PANCHROMATIC_SAMPLES;
    case lower('BQA')
        BQAfile = fullfile(dataFolder,MTL_list.PRODUCT_METADATA.FILE_NAME_BAND_QUALITY);
        LS8_BQA = geotiffread(BQAfile);
        if nargout>numOut
            varargout{1} = geotiffinfo(BQAfile);
        end
        X = unpackLandsat8BQA(LS8_BQA,datatype);
    case 'allbands'
        X_30m = GetLandsat8(dataFolder,'band1',varargin{:});
        for b=2:11
            if b == 8 %pan band is a different sized matrix
                X_15m = GetLandsat8(dataFolder,'band8',varargin{:});
                continue
            end
            X_30m = cat(3,X_30m,GetLandsat8(dataFolder,['band' num2str(b)],varargin{:}));
        end
        %create output structure
        X = struct('OLI_30m',X_30m,'OLI_15m',X_15m);
        return
    otherwise
        if strncmpi(whichVariable,'band',4)
            bandNo = whichVariable(5:end);
            band_param = 'FILE_NAME_BAND_';
            bf = strcat(band_param,bandNo);
            bandFile = fullfile(pathstr, MTL_list.PRODUCT_METADATA.(bf));
            Qcal = geotiffread(bandFile);
            
            if nargout>numOut
                varargout{1} = geotiffinfo(bandFile);
            end
            
            X = convertLandsat8DN(Qcal, MTL_list, bandNo);
        else
            error('Variable %s unrecognized',whichVariable)
        end
end
end


%if a folder, make sure last char is / or \.
function s=foldEnd(c)
if c(end) ~= '\' || c(end) ~= '/'
    g = strfind(c,'\');
    if isempty(g)
        h = strfind(c,'/');
        if ~isempty(h) && c(h(end)) ~= c(end)
            s = strcat(c,'/');
        end
    elseif c(g(end)) ~= c(end)
        s = strcat(c,'\');
    end
else
    s = c;
end
end

%extract info from filepath
function [pathstr, filetype, sceneID] = filepathInfo(filepath)
if isdir(filepath)
    filepath = foldEnd(filepath);
end


[pathstr,fileName,ext] = fileparts(filepath);
if ~isempty(fileName)
    file = strcat(fileName,ext);
    fileInfo = infoFromLandsatFilename(file);
else
    fileInfo = infoFromLandsatFilename(pathstr);
end
filetype = fileInfo.filetype;
sceneID = char(fileInfo.sceneID);
end

%find MTL filename and extract scene metadata
function [MTL_list, MTL_filename, dataFolder] = GetMTL(filepath, outputFolder)
%handle input file according to filetype

[pathstr, filetype, sceneID] = filepathInfo(filepath);

if ~isdir(filepath)
    
    
    %if tar file, untar, place contents in output folder.
    if contains(filepath,'.tar')
        [outputFolder] = extractLandsat8(filepath, outputFolder);
        %find meatdata file
        MTL_filename = fullfile(outputFolder, strcat(sceneID, '_MTL.txt'));
        %update filepath to output path
        dataFolder = foldEnd(outputFolder);
        
        
        %if band geotiff, find associated MTL metadata file
    elseif contains(filepath,'.tif')
        MTL_filename = fullfile(pathstr, strcat(sceneID, '_MTL.txt'));
        dataFolder = pathstr;
        
        %if meatdata file, continue to switch statement
    elseif contains(filepath,'_MTL.txt')
        MTL_filename = filepath;
        dataFolder = pathstr;
    end
    
    %if folder, set MTL file as the "file" the function is using
elseif strcmpi(filetype,'folder')
    MTL_filename = fullfile(filepath,strcat(sceneID,'_MTL.txt'));
    dataFolder = pathstr;
    
else
    warning('not a valid landsat file')
end

%make sure MTL file exists
assert(exist(MTL_filename, 'file') == 2, 'MTL file not found for Landsat8 scene %s',sceneID);

%Create structure of landsat scene metadata
MTL_list = MTL_parser(MTL_filename);
end