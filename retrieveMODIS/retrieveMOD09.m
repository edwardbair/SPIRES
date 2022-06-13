function  retrieveMOD09(TOKEN,folder,prefix,tile,thisYear,varargin )
% [date7dig,sortedFiles] = retrieveMOD09(folder,prefix,tile,thisYear,varargin )
% retrieves .hdf files from USGS for specified date ranges (specified by
% year, and optionally first and last days)
%
% input
%   TOKEN - may need to be generated, see getLAADSmod
%   folder - where the MODIS files live, or where they are to be placed
%   prefix - MODIS data source, e.g. MOD09GA
%   tile - location tile that we want for this processing
%   thisYear - calendar year of files desired
% optional input, in order
%   firstday, nDays - range of calendar days for which conversion and
%       transfer are wanted (default all days available for that year)
%
% output
%   date7dig - list of dates for which the .hdf files are available
%   sortedFiles - list of files

numArgs = 4;

if nargin<numArgs
    if isdeployed
        disp(['usage: ' mfilename ' folder prefix tile thisYear [, firstday nDays ]'])
    else
        disp(['usage: ' mfilename '(folder,prefix,tile,thisYear, [,firstday,nDays])'])
    end
    date7dig = [];
    sortedFiles = '';
    return
end

% parse input
p = inputParser;
default_firstDay = 1;
default_nDays = 366;
validationFcn = @(x) isnumeric(x) || ischar(x);
addRequired(p,'folder',@ischar);
addRequired(p,'prefix',@ischar);
addRequired(p,'tile',@ischar);
addRequired(p,'thisYear',validationFcn);
addOptional(p,'firstday',default_firstDay,validationFcn);
addOptional(p,'nDays',default_nDays,validationFcn);
parse(p,'folder','prefix','tile','thisYear',varargin{:})

% use local folder or Azure file share
[folder,diaryFolder] = identifyFolders(folder);
h = StartAzureDiary(mfilename,diaryFolder); %#ok<NASGU>

if ~isnumeric(thisYear)
    thisYear = str2double(thisYear);
end

% list of files
hdffiles = listhdfFiles(folder,prefix,thisYear,tile);

% process optional arguments
if ~isnumeric(p.Results.firstday)
    firstday = str2double(p.Results.firstday);
else
    firstday = p.Results.firstday;
end
if ~isnumeric(p.Results.nDays)
    nDays = str2double(p.Results.nDays);
else
    nDays = p.Results.nDays;
end
if eomday(thisYear,2)==29
    maxLastDay = 366;
else
    maxLastDay = 365;
end
lastday = min(firstday+nDays-1,maxLastDay);

disp([mfilename ' ' folder ' ' prefix ' ' tile ' '...
    num2str(thisYear)])
if nargin>numArgs
    disp(['days ' num2str(firstday) ' to ' num2str(lastday)])
end

% get the list of dates available
if ~isempty(hdffiles)
    [dateAlready,sortedFiles] = availableMODISDates(hdffiles,firstday,lastday);
else
    dateAlready = [];
    sortedFiles = '';
end
% generate list of dates needed
md1 = iso2datenum(thisYear*1000+firstday);
md2 = iso2datenum(thisYear*1000+lastday);
matdate = (md1:md2)';
possibleDate = datenum2iso(unique(min(matdate):max(matdate))',7);
disp('possible dates')
disp([possibleDate(1) possibleDate(end)])

% download missing files from USGS
if ~isequal(dateAlready,possibleDate)
    for k=1:length(possibleDate)
        thisDate = regexp(sortedFiles,...
            [prefix '\S*.A' num2str(possibleDate(k)) '\S*' tile]);
        if isempty(cell2mat(thisDate))
            downloadS = getLAADSmod(TOKEN,folder, prefix, possibleDate(k), tile );
            if downloadS.statusM~=0
                warning(['bad status return from ' mfilename ' date '...
                    num2str(possibleDate(k))])
                disp(downloadS)
            else
                disp([mfilename ': file for date ' num2str(possibleDate(k)) ' downloaded'])
            end
        else
            disp([mfilename ': file for date ' num2str(possibleDate(k)) ' already there'])
        end
    end
else
    disp(['all dates already in ' folder])
end

end