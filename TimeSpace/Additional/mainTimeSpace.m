function varargout = mainTimeSpace(snowFolder, cloudFolder, tile, thisYear, startDay, nDays, varargin)
% time-space interpolation and smoothing of multi-day cube
% [optional [outputfile,structure,original]] = mainTimeSpace(snowFolder, cloudFolder, tile, thisYear, startDay, nDays, ...)
%
% Note: uses parfor in the persistence filter and in smoothColumns
%
% Input
%   snowFolder - where the MODSCAG (or other method) RawEndmember.h5 files live
%       Each daily file is the output from JPL2H5 and contains
%       fractional endmembers for snow, vegetation, soil, shade, etc.
%   cloudFolder - where the weightCloud.h5 files live
%       Each daily file is the output from MODIScloud2h5
%   tile - MODIS tile in form hNNvNN
%   thisYear - calendar year
%   startDay - starting day of year
%   nDays - number of days in the output cubes
% Optional input, name-value pairs
%   'mask' - area name, e.g. 'sierra'
%   'cloudPersist' - maximum persistence of clouds, days (default 8)
%   'snowPersist' - minimum persistence of snow, days (default 2)
%   'dayBuffer' - for internal processing, how many days on either end to
%       avoid edge effects, default 4
%   'threshold' - value below which endmember fraction set to 0 (default 0.1)
%   'debug' - small size (typically 100-300) to test
%
% Output, .mat file with the following information
%   hdr - information about projection and referencing matrix
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
% NOTE: we're interested in snow only, so everywhere snow fraction is zero,
% other variables are set to NaN
%
% Jeff Dozier, 2015-09-08

%% syntax to run if called with too few arguments
numArgs = 6;
if nargin<numArgs
    disp(['minimum ' num2str(numArgs) '(>' num2str(nargin) ') arguments required'])
    if isdeployed
        disp('usage: mainTimeSpace snowFolder cloudFolder tile thisYear startDay nDays  ...)')
    else
        disp('usage: mainTimeSpace(snowFolder,cloudFolder,tile,thisYear,startDay,nDays, ...)')
    end
    disp('see documentation about optional arguments')
    return
end

fulltime = tic;

%% process required & optional arguments
p = inputParser;
defaultCloudPersist = 8;
defaultBuffer = 4;
defaultSnowPersist = 2;
defaultThreshold = 0.1;
defaultSmallSize = [];
validationFcn = @(x) isnumeric(x) || ischar(x);
addRequired(p,'snowFolder',@ischar)
addRequired(p,'cloudFolder',@ischar)
addRequired(p,'tile',@ischar)
addRequired(p,'thisYear',validationFcn)
addRequired(p,'startDay',validationFcn)
addRequired(p,'nDays',validationFcn)
addParameter(p,'dayBuffer',defaultBuffer,validationFcn)
addParameter(p,'cloudPersist',defaultCloudPersist,validationFcn)
addParameter(p,'snowPersist',defaultSnowPersist,validationFcn)
addParameter(p,'threshold',defaultThreshold,validationFcn)
addParameter(p,'debug',defaultSmallSize,validationFcn)
addParameter(p,'mask','',@ischar)
parse(p,snowFolder, cloudFolder, tile, thisYear, startDay, nDays, varargin{:})

%% paths to folders, and set diary folder if running deployed
% use folder or the Azure File Share?
[snowFolder,diaryFolder] = identifyFolders(snowFolder);
cloudFolder = identifyFolders(cloudFolder);
if isdeployed
    h = StartAzureDiary(mfilename,diaryFolder,...
        [tile '_year' num2str(thisYear) '_startDay' num2str(startDay,'%03d')]); %#ok<NASGU>
end
assert(exist(snowFolder,'dir')==7,'snowFolder (%s) does not exist',snowFolder)
assert(exist(cloudFolder,'dir')==7,'cloudFolder (%s) does not exist',cloudFolder)

%% char to numeric where needed
if ~isnumeric(thisYear)
    thisYear = str2double(thisYear);
end
if ~isnumeric(startDay)
    startDay = str2double(startDay);
end
if ~isnumeric(nDays)
    nDays = str2double(nDays);
end
if isnumeric(p.Results.dayBuffer)
    dayBuffer = p.Results.dayBuffer;
else
    dayBuffer = str2double(p.Results.dayBuffer);
end
if isnumeric(p.Results.cloudPersist)
    cloudPersist = p.Results.cloudPersist;
else
    cloudPersist = str2double(p.Results.cloudPersist);
end
if isnumeric(p.Results.snowPersist)
    snowPersist = p.Results.snowPersist;
else
    snowPersist = str2double(p.Results.snowPersist);
end
if isnumeric(p.Results.threshold)
    threshold = p.Results.threshold;
else
    threshold = str2double(p.Results.threshold);
end
if isempty(p.Results.debug)
    debugOption = false;
else
    debugOption = true;
    if isnumeric(p.Results.debug)
        smallSize = p.Results.debug;
    else
        smallSize = str2double(p.Results.debug);
    end
end
mask = p.Results.mask;
% MATLAB date
startDate = iso2datenum(thisYear*1000+startDay);

% warning message if missing files (doesn't abort, as a few files can be
% missing)
missingTbl = checkTimeSpaceFiles(snowFolder,cloudFolder,tile,thisYear,startDay,nDays,false);
if ~isempty(missingTbl)
    warning('some MODSCAG and/or cloud files missing')
    disp(missingTbl)
end

%% structure with raw cubes
if debugOption
    S = prepareRawCube(snowFolder,cloudFolder,tile,startDate,nDays,mask,...
        'cloudPersist',cloudPersist,'snowPersist',snowPersist,...
        'dayBuffer',dayBuffer,'threshold',threshold,'debug',smallSize);
else
    S = prepareRawCube(snowFolder,cloudFolder,tile,startDate,nDays,mask,...
        'cloudPersist',cloudPersist,'snowPersist',snowPersist,...
        'dayBuffer',dayBuffer,'threshold',threshold);
end
if isempty(S)
    warning('no cube for tile %s, startDate %s',tile,datestr(startDate));
    if nargout>0
        varargout{1} = '';
    end
    return
end
origS = S;

%% fill the NaNs, except where weight is zero
F = fillSpace(S);

%% smooth the cube along the time dimension
S = smoothTime(F);

%% make sure NaN outside mask
[endmember,attribute,~,cloudMask,imgMask,rmVariable] = identifyFields(S);
variable = cat(2,endmember,attribute);
for k=1:length(variable)
    S.(variable{k}) = setMask(S.(variable{k}),imgMask);
end
for k=1:length(rmVariable)
    S = rmfield(S,rmVariable{k});
end

%% final normalization
S.snow_fraction = snowPersistenceFilter(S.snow_fraction,snowPersist,S.threshold);
% N = normalizeEndmembers(S,endmember,S.threshold);
% for k=1:length(endmember)
%     S.(endmember{k}) = N.(endmember{k});
% end

%% set attribute values to NaN where no snow or NaN
t = isnan(S.snow_fraction) | S.snow_fraction==0;
for k=1:length(attribute)
    S.(attribute{k})(t) = NaN;
end

%% output
% truncate cubes to just the start and end days, eliminating the buffer
S.imgMask = imgMask;
S.cloudMask = cloudMask;
S = rmfield(S,'weight');
if ~isequal(S.measdate,S.allDay)
    kbegin = find(S.allDay==S.measdate(1));
    kend = find(S.allDay==S.measdate(end));
    if kbegin>=1 && kbegin<=kend && kend<=length(S.allDay);
        for k=1:length(variable)
            S.(variable{k}) = S.(variable{k})(:,:,kbegin:kend);
        end
        S.cloudMask = S.cloudMask(:,:,kbegin:kend);
    end
end
S = rmfield(S,'allDay');
file = fullfile(snowFolder,[S.prefix '_MODSCAG_A'...
    num2str(datenum2iso(S.measdate(1),7))...
    '_' num2str(datenum2iso(S.measdate(end),7))...
    '_' tile '_EndmemberCube' '.mat']);
ofile = fullfile(snowFolder,[S.prefix '_MODSCAG_A'...
    num2str(datenum2iso(S.measdate(1),7))...
    '_' num2str(datenum2iso(S.measdate(end),7))...
    '_' tile '_RawEndmemberCube' '.mat']);
save(file,'-struct','S')
save(ofile,'-struct','origS')
disp('all done')
if nargout>0
    varargout{1} = file;
    if nargout>1
        varargout{2} = S;
        if nargout>2
            varargout{3} = origS;
        end
    end
end
disp(S)
elapsedTime=toc(fulltime);
disp(['total time ' num2str(elapsedTime/3600) ' hours'])
end