function [missingTbl] =...
    checkTimeSpaceFiles(snowFolder, cloudFolder, tile, thisYear, varargin)
% [missingTbl] =checkTimeSpaceFiles(snowFolder, cloudFolder, tile, thisYear, varargin)
%checks snowFolder and cloudFolder for to make sure all files are present,
% returns list of datenums missing
%input:
%   snowFolder - where the MODSCAG (or other method) RawEndmember.h5 files live
%       Each daily file is the output from JPL2H5 and contains
%       fractional endmembers for snow, vegetation, soil, shade, etc.
%   cloudFolder - where the weightCloud.h5 files live
%       Each daily file is the output from MODIScloud2h5
%   tile - MODIS tile in form hNNvNN%   startDate - MATLAB datenum, start of period
%   thisYear - calendar year
%optional input
%   startDay - starting day of year, default 1
%   nDays - number of days in the output cubes, default rest of year
%   verbose - true or false, to print results, default true
%
% output
%   missingTbl - with dates of missing MODSCAG or cloud files
%       returns empty table if no missing files

%% syntax to run if called with too few arguments
numArgs = 4;
if nargin<numArgs
    if isdeployed
        disp(['usage: ' mfilename ' snowFolder cloudFolder tile thisYear [ startDay nDays true/false]'])
    else
        disp(['usage: ' mfilename '(snowFolder, cloudFolder, tile, thisYear, [,startDat,nDays, true/false])'])
    end
    missingTbl = [];
    return
end

%% parse input
p = inputParser;
validationFcn = @(x) isnumeric(x) || ischar(x);
v2Fcn = @(x) islogical(x) || ischar(x);
addRequired(p,'snowFolder',@ischar)
addRequired(p,'cloudFolder',@ischar)
addRequired(p,'tile',@ischar)
addRequired(p,'thisYear',validationFcn)
addOptional(p,'startDay',1,validationFcn)
addOptional(p,'nDays',366,validationFcn)
addOptional(p,'verbose',true,v2Fcn)
parse(p,snowFolder,cloudFolder,tile,thisYear,varargin{:});
snowFolder = identifyFolders(snowFolder);
cloudFolder = identifyFolders(cloudFolder);
if ~isnumeric(thisYear)
    thisYear = str2double(thisYear);
end
if isnumeric(p.Results.startDay)
    startDay = p.Results.startDay;
else
    startDay = str2double(p.Results.startDay);
end
if isnumeric(p.Results.nDays)
    nDays = p.Results.nDays;
else
    nDays = str2double(p.Results.nDays);
end
if islogical(p.Results.verbose)
    verbose = p.Results.verbose;
else
    verbose = ~strcmpi(p.Results.verbose,'false');
end

%% corresponding isodates to find the file names
startDate = iso2datenum(1000*thisYear+startDay);
startDate = floor(startDate); % eliminate parts of days, just in case
endDate = startDate+nDays-1;
% check to not go beyond end of year
[y,~,~] = datevec(endDate);
if y>thisYear;
    endDate = datenum(thisYear,12,31);
end
allDay = (startDate:endDate)';
isoDay = datenum2iso(allDay,7);

%% file names
cloudFile = cell(length(allDay),1);
modFile = cell(size(cloudFile));
modExist = false(size(allDay));
cloudExist = false(size(allDay));
assert(exist(snowFolder,'dir')==7,'snowFolder (%s) does not exist',snowFolder)
assert(exist(cloudFolder,'dir')==7,'cloudFolder (%s) does not exist',cloudFolder)
for k=1:length(allDay)
    f = dir(fullfile(snowFolder,['*_A' num2str(isoDay(k)) '*' tile '*_RawEndmembers.h5']));
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

modFileDatesMissing=allDay(~modExist);
cloudFileDatesMissing=allDay(~cloudExist);
% build table for more useful output
if any(~modExist) || any(~cloudExist)
    % table must have same set of days
    totalDays = unique(cat(1,modFileDatesMissing,cloudFileDatesMissing));
    isoDate = int32(datenum2iso(totalDays,7));
    missingMOD = cell(length(totalDays),1);
    missingCloud = cell(size(missingMOD));
    for k=1:length(totalDays)
        if any(modFileDatesMissing==totalDays(k))
            missingMOD{k} = 'missing';
        end
        if any(cloudFileDatesMissing==totalDays(k))
            missingCloud{k} = 'missing';
        end
    end
    MATLABdate = int32(round(totalDays));
    missingTbl = table(MATLABdate,isoDate,missingMOD,missingCloud);
else
    missingTbl = table;
end

%% print output if specified
if verbose
    if isempty(missingTbl)
        disp('no MODSCAG or cloud files missing')
    else
        disp('some MODSCAG and/or files missing')
        disp(missingTbl)
    end
end
end