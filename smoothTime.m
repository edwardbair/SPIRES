function [ S ] = smoothTime( U)
%smooth multi-day set of NaN-filled values
% [ smoothed ] = modscagSmoothTime( U )
%
% Input
%   U - structure with NaN-filled but unsmoothed endmembers and attributes,
%       along with confidence weights
%
% Jeff Dozier, 2016-01-16 revised to use smoothDataCube

%% process required argument
p = inputParser;
addRequired(p,'U',@isstruct)
parse(p,U)
verbose = true;

if verbose
    disp(['calling ' mfilename ' ' datestr(now)])
    disp(U)
    tstart = tic;
end

%% convert veg & rock fractions within their total
% U = convertVegRock(U,'forward');

%% process each through smoothDataCube, with the veg and rock fractions converted as above
method = 'smoothingspline';
[endmember,attribute,weight,~,imgMask] = identifyFields(U);
S = U;
weight(weight<=0 | isnan(weight)) = 0.01;
for k = 1:length(endmember)
    if verbose
        tic
    end
    % (don't time-smooth the other_fraction)
    if strcmp(endmember{k},'other_fraction')
        S.other_fraction = U.other_fraction;
        
        % note potential trouble here: the code depends on the
        % endmember list haveing snow first and vegetation before rock
    elseif strcmp(endmember{k},'rock_fraction')
        S.(endmember{k}) = 1-S.vegetation_fraction;
    elseif strcmp(endmember{k},'snow_fraction')
        X = smoothDataCube(U.(endmember{k}),weight,'mask',imgMask,'method',method);
        X(X<S.threshold) = 0;
        S.(endmember{k}) = X;
        snowMask = nansum(X,3)>0; % 2D snow mask
    else % vegetation - only smooth where there's snow
        S.(endmember{k}) = smoothDataCube(U.(endmember{k}),weight,'mask',snowMask,'method',method);
    end
    
    if verbose
        disp(['smoothed variable ' endmember{k} ' ' datestr(now)])
        toc
    end
end

for k=1:length(attribute)
    if verbose
        tic
    end
    X = U.(attribute{k});
    X(S.snow_fraction==0 | isnan(S.snow_fraction) | isnan(X)) = nanmedian(X(:));
    S.(attribute{k}) = smoothDataCube(X,weight,'mask',snowMask,'method',method);
    X = S.(attribute{k});
    X(S.snow_fraction==0 | isnan(S.snow_fraction)) = NaN;
    S.(attribute{k}) = X;
    if verbose
        disp(['smoothed variable ' attribute{k} ' ' datestr(now)])
        toc
    end
end

%% renormalize the vegetation and rock fractions
% S = convertVegRock(S,'inverse');

%% normalize endmembers
% N = normalizeEndmembers(S,endmember,S.threshold);
% for k=1:length(endmember)
%     S.(endmember{k}) = N.(endmember{k});
% end

%% done
if verbose
    elapsedTime = toc(tstart);
    disp(['total time for ' mfilename ' ' num2str(elapsedTime/3600) ' hours'])
end
end