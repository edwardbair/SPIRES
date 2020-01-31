function output_file = consolidateTimeSpace(folder,prefix,tile,filetype,varargin)
% consolidate the .mat output files from mainTimeSpace into a .mat or .h5 file
% that contains all endmembers and attributes for either the whole water year,
% the whole calendar year, calendar year and month, or arbitrary time period
%
% outputfile = consolidateTimeSpace(folder,prefix,tile,filetype, 'WY', NNNN)
% outputfile = consolidateTimeSpace(folder,prefix,tile,filetype, 'CY', NNNN)
% outputfile = consolidateTimeSpace(folder,prefix,tile,filetype, 'CYM', NNNNMM)
% outputfile = consolidateTimeSpace(folder,prefix,tile,filetype,'start',YYYYDDD,'end',YYYYDDD)
% outputfile = consolidateTimeSpace(folder,prefix,tile,filetype,'start',YYYYMMDD,'end',YYYYMMDD)
% outputfile = consolidateTimeSpace(...,'variable',variableNo)
%
% Input
%   folder - where all the output files from mainTimeSpace live (if on the
%       Azure File Share, omit the drive #)
%   prefix - typically 'MOD09GA' or 'MYD09GA'
%   tile - MODIS tile in form 'hNNvNN'
%   output filetype - either 'mat' or 'h5'
% Required input about dates
%   either 'WY' followed by NNNN, 'CY' followed by NNNN, 'start' and 'end'
%   pairs
%   (dates are ISO, can be 7- or 8-digit)
%
% Output
%   file name of output file

% endmembers have to persist for at least nPersist days above threshold
nPersist = 2;
threshold = 0.1;
output_file = '';

numArg = 6;
if nargin<numArg
    output_file = '';
    if isdeployed
        disp('if folder is on the Azure File Share, omit the drive #')
        disp('When running in deployed mode, no parentheses, no commas, no quotes.')
    end
    disp(['outputfile = ' mfilename '(folder,prefix,tile,filetype, ''WY'', NNNN)'])
    disp(['outputfile = ' mfilename '(folder,prefix,tile,filetype, ''CY'', NNNN)'])
    disp(['outputfile = ' mfilename '(folder,prefix,tile,filetype, ''CYM'', NNNNMM)'])
    disp(['outputfile = ' mfilename '(folder,prefix,tile,filetype, ''start'',YYYYDDD,''end'',YYYYDDD)'])
    disp('(dates can also be 8-digit ISO dates, YYYYMMDD)')
    return
end

totalTime = tic;

p = inputParser;
defaultWY = [];
defaultCY = [];
defaultCYM = [];
defaultStartDate = [];
defaultEndDate = [];
validationFcn = @(x) isnumeric(x) || ischar(x);
addRequired(p,'folder',@ischar)
addRequired(p,'prefix',@ischar)
addRequired(p,'tile',@ischar)
addRequired(p,'filetype',@ischar)
addParameter(p,'WY',defaultWY,validationFcn)
addParameter(p,'CY',defaultCY,validationFcn)
addParameter(p,'CYM',defaultCYM,validationFcn)
addParameter(p,'startDate',defaultStartDate,validationFcn)
addParameter(p,'endDate',defaultEndDate,validationFcn)
parse(p,folder,prefix,tile,filetype,varargin{:})
folder = p.Results.folder;
prefix = p.Results.prefix;
tile = p.Results.tile;
filetype = p.Results.filetype;

% use a local folder or the Azure File Share
[folder,diaryFolder] = identifyFolders(folder);

% which date option
waterYear = ~isempty(p.Results.WY);
calendarYear = ~isempty(p.Results.CY);
calendarYearMonth = ~isempty(p.Results.CYM);
customDates = ~(isempty(p.Results.startDate) || isempty(p.Results.endDate));
% can specify only one
howMany = [waterYear calendarYear calendarYearMonth customDates];
if nnz(howMany)~=1
    StartAzureDiary(mfilename,diaryFolder);
    error('must specify exactly one option -- ''WY'', ''CY'', ''CYM'' or ''start'',''end'' pair')
end
whichYear = find(howMany);

switch whichYear
    case 1 % water year
        if isnumeric(p.Results.WY)
            WY = p.Results.WY;
        else
            WY = str2double(p.Results.WY);
        end
        startDate = datenum([WY-1 10 1]);
        endDate = datenum([WY 9 30]);
        yr = WY;
        dstring = ['_' tile '_WY' num2str(yr)];
    case 2 % calendar year
        if isnumeric(p.Results.CY)
            CY = p.Results.CY;
        else
            CY = str2double(p.Results.CY);
        end
        startDate = datenum([CY 1 1]);
        endDate = datenum([CY 12 31]);
        yr = CY;
        dstring = ['_' tile '_CY' num2str(yr)];
    case 3 % calendar year and month
        if isnumeric(p.Results.CYM)
            cym = p.Results.CYM;
        else
            cym = str2double(p.Results.CYM);
        end
        yr = floor(cym/100);
        mon = mod(cym,100);
        startDate = datenum(yr,mon,1);
        endDate = datenum(yr,mon,eomday(yr,mon));
        mn = month(datetime(startDate,'convertfrom','datenum'),'shortname');
        monthName = char(mn);
        dstring = ['_' tile '_CY' num2str(yr) '_' monthName];
    case 4 % custom dates
        if isnumeric(p.Results.startDate)
            startDate = iso2datenum(p.Results.startDate);
            endDate = iso2datenum(p.Results.endDate);
        else
            startDate = iso2datenum(str2double(p.Results.startDate));
            endDate = iso2datenum(str2double(p.Results.endDate));
        end
        dstring = ['_' tile];
end
StartAzureDiary(mfilename,diaryFolder,dstring);

if strcmpi(filetype,'mat')
    useMAT = true;
elseif strcmpi(filetype,'h5')
    useMAT = false;
else
    error('filetype ''%s'' not recognized, must be either ''mat'' or ''h5''',filetype)
end

% list of files, which we'll then check for date range
f = dir(fullfile(folder,[prefix '*_EndmemberCube.mat']));
assert(~isempty(f),'no files that match ''%s'' in folder %s',...
    [prefix '*' tile '*_EndmemberCube.mat'], folder);
for k=1:length(f)
    possibleFiles{k,1} = fullfile(folder,f(k).name); %#ok<AGROW>
end

% date etc info from file names, then downselect to the date range
possibleTbl = infoFromMODISfilename(possibleFiles);
kfirst = find(possibleTbl{:,'matdate'}<=startDate,1,'last');
if isempty(kfirst)
    kfirst = 1;
end
klast = find(possibleTbl{:,'enddate'}>=endDate,1,'first');
if isempty(klast)
    klast = height(possibleTbl);
end

Tbl = possibleTbl(kfirst:klast,:);
% revise start/end dates based on data availability
buffer = 16;
sd = max(startDate,Tbl{1,'matdate'});
ed = min(endDate,Tbl{end,'enddate'});
if sd-startDate>buffer
    startDate = sd;
end
if endDate-ed>buffer
    endDate = ed;
end
allDays = (startDate:endDate)';

disp('file names and their information')
disp(['startDate ' datestr(startDate) ', endDate ' datestr(endDate)])
disp(Tbl)

% put together the cubes for each endmember in sequence
% read the first file, and specify which variables in this sequence
m1 = matfile(char(Tbl{1,'file'}));
[variable,snowAtt,DS] = specifyVariables(m1);
startNo = 1;
finNo = length(variable);

for e=startNo:finNo
    tic
    disp(['processing ' variable{e}])
    % first file has been accessed, do the others to process and
    % derive the output file name
    [output_file, dStruct] = nameFile(useMAT,tile,DS,folder,prefix,dstring,...
        startDate,endDate);
    X = m1.(variable{e});
    mdX = m1.measdate;
    mdX = (mdX(:))';
    kfirst = find(mdX>=startDate,1,'first');
    klast = find(mdX<=endDate,1,'last');
    assert(~isempty(kfirst) && ~isempty(klast),...
        'something wrong, kfirst or klast empty')
    disp(['first date ' datestr(mdX(kfirst))])
    consolX = X(:,:,kfirst:klast);
    consolD = mdX(kfirst:klast);
    for k=2:height(Tbl) % next file in the date sequence
        disp(['     processing ' char(Tbl{k,'file'})])
        mk = matfile(char(Tbl{k,'file'}));
        Y = mk.(variable{e});
        mdY = mk.measdate;
        mdY = (mdY(:))';
        kfirst = find(mdY>=startDate,1,'first');
        klast = find(mdY<=endDate,1,'last');
        consolX = cat(3,consolX,Y(:,:,kfirst:klast));
        consolD = cat(2,consolD,mdY(kfirst:klast));
    end
    disp(['last date ' datestr(mdY(klast))])
    toc
    % end members have to persist for at lease nPersist days
    isMember = true;
    % don't consider "attributes," which are distinct from "endmembers"
    for a=1:length(snowAtt)
        if strcmp(variable{e},snowAtt{a})
            isMember = false;
            break;
        end
    end
    % fill if missing days
    if ~isequal(consolD,allDays)
        consolX = fillGaps(consolX,consolD,allDays);
    end
    % persistence
    if isMember
        consolX = snowPersistenceFilter(consolX,nPersist,threshold);
    end
    cube2file(output_file,dStruct,allDays,variable{e},consolX);
    disp([output_file ' written with endmember ' variable{e}])
    toc
end
disp(['total elapsed time ' num2str(toc(totalTime)/3600) ' hours'])
end

function Z = fillGaps(Y,origX,totalX)
% transpose to put columns in contiguous memory
N = size(Y);
Y = reshape(Y,N(1)*N(2),N(3))';
Z = nan(length(totalX),size(Y,2));
parfor c=1:size(Y,2)
    if ~all(isnan(Y(:,c)))
        y = Y(:,c);
        t = ~isnan(y);
        % extrapolate to the ends with nearest neighbors
        m = find(t,1,'first');
        if m~=1
            y(1:m-1) = y(m);
        end
        m = find(t,1,'last')
        if m<length(y)
            y(m+1:end) = y(m);
        end
        t = ~isnan(y);
        Z(:,c) = interp1(origX(t),y(t),totalX,'linear'); %#ok<PFBNS>
    end
end
Z = reshape(Z',N(1),N(2),length(totalX));
end

function [file,dStruct] = nameFile(useMAT,tile,DS,folder,prefix,dstring,...
    startDate,endDate)
% name the output file
% Output
%   file - full path to file
%   dStruct - structure with information about projection, etc

persistent already filename DShold

if isempty(already)
    already = true;
    
    d1 = num2str(datenum2iso(startDate,7));
    d2 = num2str(datenum2iso(endDate,7));
    
    if useMAT
        suffix = '.mat'; %#ok<NASGU>
    else
        suffix = '.h5'; %#ok<NASGU>
    end
    
    % information common to all datasets, not used in this
    % program but want in the output file
    [R,ProjectionStructure] = sinusoidProjMODtile(tile);
    DS.RefMatrix = R.RefMatrix_500m;
    DS.ProjectionStructure = ProjectionStructure;
    DShold = DS;
    
    if useMAT
        suffix = '.mat';
    else
        suffix = '.h5';
    end
    filename = [prefix '_TimeSpaceCube_A' d1 '_' d2 dstring  suffix];
end

% these values returned even if the routine has already been called
dStruct = DShold;
file = fullfile(folder,filename);
end

function cube2file(filename,dStruct,matdate,member,Value)
% write this endmember into the file, which is created if it doesn't
% already exist
persistent already

% which kind of file
if ~isempty(strfind(filename,'.mat'))
    useMAT = true;
elseif ~isempty(strfind(filename,'.h5'))
    useMAT = false;
else
    error('filename %s, unrecognized extension',filename)
end

% if this member is deltavis, range should be 0-1 not 0-100
if strcmp(member,'deltavis') && max(Value(:))>1
    Value = Value/100;
end

% convert to scaled integer
D = dStruct;
X.(member) = float2integer(Value,D.(member).divisor,0,D.(member).dataType,0,D.(member).maxVal);
if isfield(D.(member),'FillValue')
    if D.(member).FillValue~=0
        X.(member)(isnan(Value)) = D.(member).FillValue;
    end
end

if useMAT % .mat file
    if isempty(already)
        already = true;
        if exist(filename,'file')
            delete(filename)
        end
        RefMatrix = dStruct.RefMatrix; %#ok<NASGU>
        ProjectionStructure = dStruct.ProjectionStructure; %#ok<NASGU>
        dStruct = rmfield(dStruct,'RefMatrix');
        dstruct = rmfield(dStruct,'ProjectionStructure'); %#ok<NASGU>
        ISOdates = datenum2iso(matdate,7); %#ok<NASGU>
        MATLABdates = matdate; %#ok<NASGU>
        save(filename,'RefMatrix','ProjectionStructure','dStruct',...
            'MATLABdates','ISOdates','-v7.3')
    end
    
    % save cube to file
    save(filename,'-struct','X','-append')
    
else % .h5 file
    if isempty(already)
        if exist(filename,'file')
            delete(filename)
        end
    end
    group = '/Grid/MODIS_GRID_500m'; % MODSCAG data at 500 m resolution
    arraySize = size(Value);
    chunkSize = [arraySize(1) arraySize(2) 1];
    deflateLevel = 9;
    % write data to file
    if isfield(D.(member),'FillValue')
        h5create(filename,[group '/' member],arraySize,...
            'Deflate',deflateLevel,...
            'ChunkSize',chunkSize,...
            'DataType',D.(member).dataType,...
            'FillValue',D.(member).FillValue)
    else
        h5create(filename,[group '/' member],arraySize,...
            'Deflate',deflateLevel,...
            'ChunkSize',chunkSize,...
            'DataType',D.(member).dataType)
    end
    h5write(filename,[group '/' member],X.(member))
    h5writeatt(filename,[group '/' member],'divisor',D.(member).divisor)
    if isfield(D.(member),'units')
        h5writeatt(filename,[group '/' member],'units',D.(member).units)
    end
    
    % initial values
    if isempty(already)
        already = true;
        % referencing matrix and projection information
        h5writeProjection(filename,'/Grid',dStruct.ProjectionStructure)
        h5writeatt(filename,group,'ReferencingMatrix',dStruct.RefMatrix)
        % write dates
        ISOdates = datenum2iso(matdate,7);
        MATLABdates = matdate;
        h5writeatt(filename,'/','MATLABdates',MATLABdates);
        h5writeatt(filename,'/','ISOdates',ISOdates);
    end
end
end

function [variable,snowAtt,DS] = specifyVariables(m1)
possible = {'snow_fraction','vegetation_fraction','rock_fraction',...
    'other_fraction','deltavis','grain_size','drfs_grnsz'};
snowAtt = {'grain_size','drfs_grnsz','deltavis'}; % (not really 'endmembers')
sizeUnit = 'mu_m';
mfield = fieldnames(m1);
% structure with information on datatypes etc.
% possible variables and their attributes
dataType = 'int16';
maxVals = [ones(1,4) 0.6 ones(1,2)*3200];
divisor = [ones(1,5)*10000 ones(1,2)*10];
v = 0;
for a=1:length(possible)
    list = strcmp(possible{a},mfield);
    if ~isempty(find(list,1))
        v = v+1;
        variable{v} = possible{a}; %#ok<AGROW>
        DS.(variable{v}).dataType = dataType;
        DS.(variable{v}).maxVal = maxVals(a);
        DS.(variable{v}).divisor = divisor(a);
        DS.(variable{v}).FillValue = intmin(dataType);
        if strcmp('grain_size',variable{v}) || strcmp('drfs_grnsz',variable{v})
            DS.(variable{v}).units = sizeUnit;
        end
    end
end
end