function [ group ] = findMODISh5group( filename,pattern )
%find group in MODIS h5 file that matches a pattern
% [ group ] = findMODISh5group( filename,pattern )
%   (this function needed because of historical changes in naming
%   convention)
%
% Input
%   filename - MODIS-related .h5 filename
%   pattern - to match, e.g. '500m'
%   (assumed to be in the '/Grid' part of the file)
%
% Output
%   group - name of group, e.g. '/Grid/MODIS_Grid_500m'

info = h5info(filename,'/Grid');
group = '';
for k=1:length(info.Groups)
    if strfind(info.Groups(k).Name,pattern)
        group = info.Groups(k).Name;
        break;
    end
end
assert(~isempty(group),'%s group not found in file %s',pattern,filename)

end