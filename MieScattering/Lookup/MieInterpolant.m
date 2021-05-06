function [ F,inputVariable,outputVariable ] = MieInterpolant( C, M, thisVariable )
% [ F,inputVariable,outputVariable ] = MieInterpolant( C, M, thisVariable )
%create interpolating function based on results from prepMieInterpolant
%
%Input
%   C - structure with inputs as N-dimensional grids
%   M - Mie structure with outputs as N-D grids
%
%Output
%   F's - cell vector of gridded interpolant objects for these inputs
%   variable - names of predicted variable in order of the F cell vector

%Variable names (in the structures, they're N-D objects)
if strcmpi(thisVariable,'wetsnow')
    inputVariable = {'waterFraction','wavelength','radius'};
    wetSnow = true;
else
    inputVariable = {'wavelength','radius'};
    wetSnow = false;
end
outputVariable = {'omega','g','Qext'}; % don't need the other Q's but could add

method = 'makima';
extrapolation = 'nearest';

holdF = cell(length(outputVariable),1);
for k=1:length(outputVariable)
    % for values close to 1.0, i.e. single-scattering albedo, use log(1 minus)
    if strcmp(outputVariable{k},'omega')
        y = log(1-M.(outputVariable{k}));
    else
        y = M.(outputVariable{k});
    end
    if wetSnow
        F = griddedInterpolant(C.waterF,log(C.wavelength),...
            sqrt(C.radius),y,method,extrapolation);
    else
        F = griddedInterpolant(log(C.wavelength),sqrt(C.radius),...
            y,method,extrapolation);
    end
    holdF{k} = F;
end

F = holdF;

end