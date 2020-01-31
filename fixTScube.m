function output_file = fixTScube(folder,prefix,tile,WY)
% output_file = fixTScube(folder,prefix,tile,WY)
% fix the TS cubes to make sure they contain the whole year
% (also need to fix consolidateTimeSpace to do this)
%
% Input
%   folder - where all the current files live (if on the
%       Azure File Share, omit the drive #)
%   prefix - typically 'MOD09GA' or 'MYD09GA'
%   tile - MODIS tile in form 'hNNvNN'
%   WY - water year
%
% Output
%   file name of output file

numArg = 4;
if nargin<numArg
    output_file = '';
    if isdeployed
        disp('if folder is on the Azure File Share, omit the drive #')
        disp('When running in deployed mode, no parentheses, no commas, no quotes.')
    end
    disp(['outputfile = ' mfilename '(folder,prefix,tile,WY)'])
    return
end

totalTime = tic;

p = inputParser;
validationFcn = @(x) isnumeric(x) || ischar(x);
addRequired(p,'folder',@ischar)
addRequired(p,'prefix',@ischar)
addRequired(p,'tile',@ischar)
addRequired(p,'WY',validationFcn)
parse(p,folder,prefix,tile,WY)
folder = p.Results.folder;
prefix = p.Results.prefix;
tile = p.Results.tile;
if isnumeric(p.Results.WY)
    WY = p.Results.WY;
else
    WY = str2double(p.Results.WY);
end
% use a local folder or the Azure File Share
[folder,diaryFolder] = identifyFolders(folder);
StartAzureDiary(mfilename,diaryFolder,['WY' num2str(WY)]);

% file name
f = dir(fullfile(folder,[prefix '*TimeSpaceCube*' tile '*WY' num2str(WY) '*.h5']));
if length(f)~=1
    disp(['# files = ' num2str(length(f))])
    disp(f)
    error('should be just one file')
end

% file information
fname = fullfile(folder,f.name);
I = h5info(fname);
T = infoFromMODISfilename(fname);
allDays = (T.matdate:T.enddate)';
disp(['start date ' num2str(datenum2iso(T.matdate,7)) ' enddate '...
    num2str(datenum2iso(T.enddate,7))])
if length(I.Attributes(1).Value)==length(allDays)
    disp(['file ' f.name ' has all days, stopping'])
    output_file = fname;
    return
else
    fileDays = I.Attributes(1).Value;
    if length(fileDays)>length(allDays)
        disp(['file ' f.name ' has ' num2str(length(fileDays))...
            ' days, so will truncate to ' num2str(length(allDays))])
    else
        disp(['file ' f.name ' has ' num2str(length(fileDays))...
            ' days, so will interpolate to ' num2str(length(allDays))])
    end
    output_file = fullfile(folder,['new' f.name]);
end

mfield = {'snow_fraction','grain_size','deltavis'};
[variable,~,DS] = specifyVariables(mfield);
for k=1:length(variable)
    tic
    [X,dates,hdr] = GetEndmember(fname,variable{k});
    % deal with non-unique dates (error in consolidateTimeSpace)
    if ~isequal(dates,unique(dates))
        [dates,ia] = unique(dates);
        X = X(:,:,ia);
    end
    allTogether = isequal(dates,allDays);
    DS.ProjectionStructure = hdr.ProjectionStructure;
    DS.RefMatrix = hdr.RefMatrix;
    if allTogether
        Z = X;
    else
        % transpose to put columns in contiguous memory
        N = size(X);
        X = reshape(X,N(1)*N(2),N(3))';
        x = double(dates);
        allNaN = false(1,size(X,2));
        parfor c=1:size(X,2)
            allNaN(c) = all(isnan(X(:,c)));
        end
        whichCols = find(~allNaN);
        newX = X(:,whichCols);
        newZ = zeros(length(allDays),size(newX,2));
        sumX = nansum(newX,1);
        parfor c=1:length(whichCols);
            if sumX(c)>0 % all zeros - leave at zero
                y = double(newX(:,c));
                t = ~isnan(y);
                try
                    newZ(:,c) = interp1(x(t),y(t),allDays,'nearest');
                catch
                    newZ(:,c) = nan(length(allDays),1);
                end
            end
        end
        Z = nan(length(allDays),size(X,2));
        Z(:,whichCols) = newZ;
        Z = reshape(Z',N(1),N(2),length(allDays));
        if strcmp(variable{k},'snow_fraction')
            noSnow = Z<=0 | isnan(Z);
        else
            Z(noSnow) = NaN;
        end
    end
    cube2file(output_file,DS,allDays,variable{k},Z);
    disp([output_file ' written with variable ' variable{k}])
    toc
end
disp(['total elapsed time ' num2str(toc(totalTime)/3600) ' hours'])
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
    if D.(member).FillValue~=0 && min(X.(member)(:))~=D.(member).FillValue
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

function [variable,snowAtt,DS] = specifyVariables(mfield)
possible = {'snow_fraction','vegetation_fraction','rock_fraction',...
    'other_fraction','deltavis','grain_size','drfs_grnsz'};
snowAtt = {'grain_size','drfs_grnsz','deltavis'}; % (not really 'endmembers')
sizeUnit = 'mu_m';
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