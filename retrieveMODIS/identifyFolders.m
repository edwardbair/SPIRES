function [ folderName,varargout ] = identifyFolders(folder )
% [ folderName [,diaryFolder] ] = identifyFolders(folder )
%identify input/output/diary folders on Windows PC, Windows Azure, Mac, or
%Unix systems
%
% Input
%   folder - full path name of input or output folder, but omit path if on
%       Azure
% Output
%   full path to folder
%   if second output argument specified, full path to diary folder

% add path if on Azure file share
if ismac || isunix
    % put hooks here for shared file system on Unix
    folderName = folder;
    az = false;
elseif ispc
    if isempty(strfind(folder,':\')) && isempty(strfind(folder,'/'))
        drive = AzureFileShare();
        folderName = fullfile(drive,folder);
        az = true;
    else
        az = false;
        folderName = folder;
    end
else
    error('neither Mac, Unix, nor Windows?')
end

% make sure folder exists
assert(exist(folderName,'file')==7,'folder %s does not exist',folderName) 

% second folder for diary
if nargout>1
    if az
        diaryFolder = fullfile(drive,'output');
        if exist(diaryFolder,'file')~=7
            mkdir(diaryFolder);
        end
    else
        diaryFolder = folder;
    end
    varargout{1} = diaryFolder;
end

end