function [ F ] = fillSpace( S)
%examine multi-day set of end member images (typically 8-32 days), interpolate
%across NaNs, and normalize
% [ F ] = fillSpace( S)
%
% Input
%   S - structure containing values from prepareRawCube, along with
%   information about dates and threshold
%       3D arrays of endmembers (which sum to 1) and attributes
%
% Output, same structure, but with NaNs filled in
%
% Jeff Dozier, 2016-01-16 revised to use new infillDataCube routine

%% process required argument
p = inputParser;
addRequired(p,'S',@isstruct)

verbose = true;
disp(['calling ' mfilename ' ' datestr(now)])
disp(S)
fillStart = tic;

%% identify the endmembers and attributes
[endmember,attribute,weight,cloudMask,imgMask] = identifyFields(S);
variable = cat(2,endmember,attribute);
for k=1:length(variable)
    S.(variable{k})(cloudMask) = NaN;
end

%% for vegetation and rock, specify the fraction of that total
% S = convertVegRock(S,'forward');

%% interpolate the NaNs
variable = cat(2,endmember,attribute);
F = S;
% snow fraction first
if verbose
    memStart=tic;
end
[F.snow_fraction] = infillDataCube(S.snow_fraction,imgMask);
W = ones(size(weight));
W(isnan(S.snow_fraction)) = 0;
F.weight = weight.*W;
F.cloudMask = cloudMask;
F.imgMask = imgMask;
if verbose
    disp(['filled variable snow_fraction ' datestr(now)])
    toc(memStart)
end
% other endmembers
for k = 1:length(endmember)
    if verbose
        memStart=tic;
    end
    if ~strcmp(endmember{k},'snow_fraction') % snow_fraction already done
        if strcmp(endmember{k},'other_fraction') % treat other_fraction differently as so few
            X = S.other_fraction;
            X(isnan(X)) = 0;
            F.other_fraction = setMask(X,F.imgMask);
        elseif strcmp(endmember{k},'rock_fraction')
            F.(endmember{k}) = 1-F.vegetation_fraction;
        else
            F.(endmember{k}) = infillDataCube(S.(endmember{k}),F.imgMask);
        end
    end
    % 2D adaptive Wiener filter, but not on other_fraction
    if ~strcmp(F.(endmember{k}),'other_fraction')
        X = F.(endmember{k});
        for d=1:size(X,3)
            X(:,:,d) = wiener2(X(:,:,d));
        end
        X(isnan(X)) = 0;
        F.(endmember{k}) = X; 
    end
    if verbose
        disp(['filled and filtered variable ' endmember{k} ' ' datestr(now)])
        toc(memStart)
    end
end
% variables associated with snow
for k=1:length(attribute)
    if verbose
        memStart=tic;
    end
    if any(F.snow_fraction(:)>0)
        % set values away from snow to median, for interpolation
        X = S.(attribute{k});
        X(F.snow_fraction==0 | isnan(F.snow_fraction)) =  nanmedian(X(:));
        F.(attribute{k}) = infillDataCube(X,F.imgMask);
    else
        F.(attribute{k}) = nan(size(F.snow_fraction),'like',F.snow_fraction);
    end
    if verbose
        disp(['filled variable ' attribute{k} ' ' datestr(now)])
        toc(memStart)
    end
end

%% renormalize the vegetation and rock fractions
% F = convertVegRock(F,'inverse');

%% normalize endmembers to sum to 1.0 and make sure attributes are NaN
%  where there's no snow
% N = normalizeEndmembers(F,endmember,F.threshold);
% for k=1:length(endmember)
%     F.(endmember{k}) = N.(endmember{k});
% end
for k=1:length(attribute)
    X = F.(attribute{k});
    X(F.snow_fraction==0 | isnan(F.snow_fraction)) = NaN;
    F.(attribute{k}) = X;
end
%% reapply the image mask, just to make sure
for k=1:length(variable)
    F.(variable{k}) = setMask(F.(variable{k}),F.imgMask);
end

if verbose
    elapsedTime = toc(fillStart);
    disp(['total time for ' mfilename ' ' num2str(elapsedTime/3600) ' hours'])
end
end