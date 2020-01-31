function [ F,variable ] = MieInterpolant( folder, template )
% [ F,variable ] = MieInterpolant( folder, template )
%create interpolating function based on master table of results from
%prepMieInterpolant
%
%Input
%   folder - where files are stored
%   template - regular expression that characterizes the files to process
%
%Output
%   F's - cell vector of gridded interpolant objects for these files
%   variable - names of predicted variable in order of the F cell vector

p = inputParser;
addRequired(p,'folder',@ischar)
addRequired(p,'template',@ischar)
parse(p,folder,template);

%list of files
f = dir(fullfile(p.Results.folder,p.Results.template));
if isempty(f)
    error('no files found that match template %s',p.Results.template)
else
    disp([num2str(length(f)) ' files found matching template ' p.Results.template])
end

% build master table
for k=1:length(f)
    m = matfile(fullfile(p.Results.folder,f(k).name));
    if k==1
        masterTbl = m.T;
    else
        masterTbl = [masterTbl; m.T]; %#ok<AGROW>
    end
end

%sort order for building interpolant
T = unique(masterTbl);

%Variable names, and turn vectors in N-D objects
variable = fieldnames(T);
wetSnow = contains(template,'wetSnow','IgnoreCase',true);
if wetSnow
    variable = variable(4:9);
    [V1,V2,V3] = ndgrid(unique(T.waterFraction),sqrt(unique(T.radius)),...
        log(unique(T.wavelength)));
else
    variable = variable(3:8);
    [V1,V2] = ndgrid(sqrt(unique(T.radius)),log(unique(T.wavelength)));
end

method = 'linear';
extrapolation = 'none';

holdF = cell(length(variable),1);
for k=1:length(variable)
    % for values close to 1.0, i.e. single-scattering albedo, use 1 minus
    if strcmp(variable{k},'omega') &&...
            (contains(template,'ice','IgnoreCase',true) ||...
            contains(template,'water','IgnoreCase',true) ||...
            contains(template,'snow','IgnoreCase',true))
        y = 1-T.(variable{k});
    else
        y = T.(variable{k});
    end
    if wetSnow
        F = griddedInterpolant(V1,V2,V3,reshape(y,size(V1)),method,extrapolation);
    else
        F = griddedInterpolant(V1,V2,reshape(y,size(V1)),method,extrapolation);
    end
    holdF{k} = F;
end

F = holdF;

end