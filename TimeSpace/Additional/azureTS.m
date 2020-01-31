function azureTS( tile, thisYear, startDay, nDays)
% azureTS( tile, thisYear, startDay, nDays)
%wrapper for time-space on Azure, for less complicated load
%
% note - uses parfor in a couple of the internal routines
% also requires the data file BoundaryShapes.mat
%
numarg = 4;
if isdeployed && nargin~=numarg
    disp(['usage: ' mfilename ' tile thisYear startDay nDays'])
    return
end
p = inputParser;
validationFcn = @(x) isnumeric(x) || ischar(x);
addRequired(p,'tile',@ischar)
addRequired(p,'thisYear',validationFcn)
addRequired(p,'startDay',validationFcn)
addRequired(p,'nDays',validationFcn)
parse(p,tile,thisYear,startDay,nDays)
% variables
if ~isnumeric(p.Results.thisYear)
    thisYear = str2double(p.Results.thisYear);
end
if ~isnumeric(p.Results.startDay)
    startDay = str2double(p.Results.startDay);
end
if ~isnumeric(p.Results.nDays)
    nDays = str2double(p.Results.nDays);
end
% mask depending on tile name
if strfind(tile,'h08')
    mask = 'sierra';
elseif strfind(tile,'h23')
    mask = 'hindukush';
else
    mask = '';
end

% folders
snowFolder = ['modscag\' tile];
cloudFolder = ['modisCloud\' tile];
if isempty(mask)
    mainTimeSpace(snowFolder, cloudFolder, tile, thisYear, startDay, nDays);
else
    mainTimeSpace(snowFolder, cloudFolder, tile, thisYear, startDay, nDays, 'mask',mask);
end

end
