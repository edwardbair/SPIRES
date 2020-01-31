function [ modS ] = prepareRawCube( snowFolder, cloudFolder, tile, startDate, nDays, mask, varargin )
%prepare multi-day set of snow endmember images (typically 8-32 days) for
%furtherprocessing
%
% [ modS ] = prepareRawCube( snowFolder,cloudFolder,tile,startDate,nDays [,'dayBuffer',value,'threshold',value])
% Input
%   snowFolder - where the MODSCAG (or other method) RawEndmember.h5 files live
%       Each daily file is the output from JPL2H5 and contains
%       fractional endmembers for snow, vegetation, soil, shade, etc.
%   cloudFolder - where the weightCloud.h5 files live
%       Each daily file is the output from MODIScloud2h5
%   tile - MODIS tile in form hNNvNN
%   startDate - MATLAB datenum, start of period
%   nDays - number of days in the output cubes
%   mask - name of mask to apply, '' if none
% Optional input, name-value pairs
%   'cloudPersist' - maximum persistence of clouds, days (default 7)
%   'snowPersist' - minimum persistence of snow, days (default 3)
%   'dayBuffer' - for internal processing, how many days on either end to
%       avoid edge effects, default 4
%   'threshold', specifying fraction below which set to zero
%       number between 0 and 1, typically not too different from zero,
%       default 0.1
%   'debug', small subset for debugging, following by linear size of square
%       region (generally in range 100 to 300)
%
% Output, modS
%   hdr - information about projection and referencing matrix
%   allDay - dates covered, in datenum format, including buffer
%   measdate - dates covered, in datenum format, excluding buffer
%   threshold - scalar, fraction below which set to zero
%   3D arrays, which could include
%       snow fraction
%       vegetation fraction
%       rock fraction
%       other fraction
%       grain size, um
%       drfs grain size, um
%       deltavis, due to dust
%       weight - confidence weights for further processing
%       cloudMask - true where cloudes
%   (missing dates within the sequence are NaN)
%
% Jeff Dozier, 2015-08-11

%% process required & optional arguments
p = inputParser;
defaultThreshold = 0.1;
defaultCloudPersist = 7;
defaultSnowPersist = 3;
defaultBuffer = 4;
defaultSmallSize = [];
addRequired(p,'snowFolder',@ischar)
addRequired(p,'cloudFolder',@ischar)
addRequired(p,'tile',@ischar)
addRequired(p,'startDate',@isnumeric)
addRequired(p,'nDays',@isnumeric)
addRequired(p,'mask',@ischar)
addParameter(p,'threshold',defaultThreshold,@isnumeric)
addParameter(p,'dayBuffer',defaultBuffer,@isnumeric)
addParameter(p,'cloudPersist',defaultCloudPersist,@isnumeric)
addParameter(p,'snowPersist',defaultSnowPersist,@isnumeric)
addParameter(p,'debug',defaultSmallSize,@isnumeric)
parse(p,snowFolder, cloudFolder, tile, startDate, nDays, mask, varargin{:})
threshold = p.Results.threshold;
dayBuffer = p.Results.dayBuffer;
if isempty(p.Results.debug)
    debugOption = false;
else
    debugOption = true;
    smallSize = p.Results.debug;
end

tic

%% corresponding isodates to find the file names
startDate = floor(startDate); % eliminate parts of days, just in case
endDate = startDate+nDays-1;
bufferStart = startDate-dayBuffer;
bufferEnd = endDate+dayBuffer;
allDay = (bufferStart:bufferEnd)';
isoDay = datenum2iso(allDay,7);

% file names
cloudFile = cell(length(allDay),1);
modFile = cell(size(cloudFile));
modExist = false(size(allDay));
cloudExist = false(size(allDay));
for k=1:length(allDay)
    f = dir(fullfile(snowFolder,['*_A' num2str(isoDay(k)) '*' tile '*RawEndmembers.h5']));
    if ~isempty(f)
        modFile{k} = fullfile(snowFolder,f.name);
        modExist(k) = true;
    end
    f = dir(fullfile(cloudFolder,['*_A' num2str(isoDay(k)) '*' tile '*_weightCloud.h5']));
    if ~isempty(f)
        cloudFile{k} = fullfile(cloudFolder,f.name);
        cloudExist(k) = true;
    end
end

% table of files present, ignore unless both cloud and modscag exist
fT = table(allDay,modExist,cloudExist,modFile,cloudFile);
t = ~(fT.modExist & fT.cloudExist);
fT(t,:) = [];
assert(height(fT)>1,...
    'not enough files for %s between %s and %s',...
    tile,datestr(startDate),datestr(endDate));
% make sure all the same prefix
ftS = infoFromMODISfilename(fT.modFile);
ftC = infoFromMODISfilename(fT.cloudFile);
assert(isequal(ftS.prefix,ftC.prefix),...
    'check data: prefixes in snow and clouds files not same')
for k=2:height(ftS)
    assert(strcmp(ftS{1,'prefix'},ftS{k,'prefix'}),...
        'check input: prefixes of files not the same')
end
prefix = ftS{1,'prefix'};
% truncate at ends
allDay = (fT{1,'allDay'}:fT{end,'allDay'})'; %#ok<BDSCA>
startDate = max(startDate,allDay(1));
endDate = min(endDate,allDay(end));
measdate = (startDate:endDate)';
if length(measdate)<2
    warning('not enough files for %s between %s and %s',...
        tile,datestr(startDate),datestr(endDate))
    modS = [];
    return
end

%% read the cloud files and build the cloud and weight cubes
[ cloudMask, weight, dateval, hdr] = GetCloudWeight(fT.cloudFile);
% are there missing dates?
missingDays = ~isequal(allDay,dateval);

% allocate the data cubes (with NaNs so any missing values are NaNs)
endmember = {'snow_fraction','vegetation_fraction','rock_fraction',...
    'other_fraction'};
attribute = {'grain_size','drfs_grnsz','deltavis'};
dsName = cat(2,endmember,attribute);
if missingDays
    for k=1:length(dsName)
        rawCube.(dsName{k}) = NaN(size(weight,1),size(weight,2),...
            length(allDay),'like',weight);
    end
    oldWeight = weight;
    weight = nan(size(rawCube.(dsName{1})));
    oldCloud = cloudMask;
    cloudMask = false(size(rawCube.(dsName{1})));
else
    for k=1:length(dsName)
        rawCube.(dsName{k}) = NaN(size(weight),'like',weight);
    end
end

% check to make sure we got the right number of days
assert(isequal(fT.allDay,dateval),...
    'cloud files don''t match table of date values')

%% build the raw cubes, by date
for d=1:height(fT)
    variableList = MODSCAGvariables(fT{d,'modFile'});
    if missingDays
        thisDay = find(allDay==fT{d,'allDay'});
        if ~isempty(thisDay)
            for k=1:length(dsName)
                list = strcmp(dsName{k},variableList);
                if nnz(list)>0
                    rawCube.(dsName{k})(:,:,thisDay) =...
                        GetEndmember(fT{d,'modFile'},dsName{k});
                end
            end
            weight(:,:,thisDay) = oldWeight(:,:,d);
            cloudMask(:,:,thisDay) = oldCloud(:,:,d);
        end
    else
        for k=1:length(dsName)
            list = strcmp(dsName{k},variableList);
            if nnz(list)>0
                rawCube.(dsName{k})(:,:,d) =...
                    GetEndmember(fT{d,'modFile'},dsName{k});
            end
        end
    end
end

%% if a mask is specified, set everything outside it to NaN, zero, or false
if ~isempty(p.Results.mask)
    [~,ProjectionStructure,RasterReference] = sinusoidProjMODtile(tile);
    imgMask = createBoundaryMask('area',p.Results.mask,'proj',...
        ProjectionStructure,...
        'rasterref',RasterReference.RasterReference_500m);
    cloudMask = setMask(cloudMask,imgMask);
    weight = setMask(weight,imgMask,0);
    for k=1:length(dsName)
        rawCube.(dsName{k}) = setMask(rawCube.(dsName{k}),imgMask);
    end
end

%% check clouds and snow for persistence
cloudMask = cloudPersistenceFilter(cloudMask,p.Results.cloudPersist);
rawCube.snow_fraction = snowPersistenceFilter(rawCube.snow_fraction,...
    p.Results.snowPersist,threshold);

%% for clouds or weights of zero, set values to NaN, and set missing weights to 0.5
weight(isnan(weight)) = 0.5;
zw = weight==0;
% 1st percentile grain size is unlikely, so set to NaN
if nnz(~isnan(rawCube.snow_fraction))>0 && isfield(rawCube,'grain_size')
    rawCube.grain_size(rawCube.snow_fraction==0) = NaN;
    t = isnan(rawCube.grain_size);
    minG = ceil(prctile(rawCube.grain_size(~t),1)*1000)/1000;
    pCloud = rawCube.grain_size<=minG;
    rawCube.grain_size(pCloud) = NaN;
end
for k=1:length(dsName)
    rawCube.(dsName{k})(zw|cloudMask) = NaN;
end
raw = normalizeEndmembers(rawCube,endmember,threshold);
for k=1:length(attribute)
    raw.(attribute{k}) = rawCube.(attribute{k});
end

%% output
modS.hdr = hdr;
modS.threshold = threshold;
modS.tile = tile;
modS.prefix = char(prefix);
modS.allDay = allDay;
modS.measdate = measdate;
modS.weight = weight;
modS.cloudMask = cloudMask;
modS.cloudPersist = p.Results.cloudPersist;
modS.snowPersist = p.Results.snowPersist;
if exist('imgMask','var')==1
    modS.imgMask = imgMask;
end
% just the raw values if they were all NaN
for k=1:length(dsName)
    if nnz(~isnan(raw.(dsName{k})))>0
        modS.(dsName{k}) = raw.(dsName{k});
    end
end

%% debug option - reduce size of cubes around area where max snow
if debugOption
    sumSnow = nansum(modS.snow_fraction,3);
    [r,c] = find(sumSnow==max(sumSnow(:)));
    fr = round(max(1,r-smallSize/2));
    lr = round(min(size(weight,1),r+smallSize/2));
    fc = round(max(1,c-smallSize/2));
    lc = round(min(size(weight,2),c+smallSize/2));
    modS.weight = weight(fr:lr,fc:lc,:);
    modS.cloudMask = cloudMask(fr:lr,fc:lc,:);
    if exist('imgMask','var')==1
        modS.imgMask=imgMask(fr:lr,fc:lc,:);
    end
    for k=1:length(dsName)
        if nnz(~isnan(raw.(dsName{k})))>0
            modS.(dsName{k}) = raw.(dsName{k})(fr:lr,fc:lc,:);
        end
    end
end
%% that's all
disp('input cubes read')
toc
end