function [ convS ] = convertVegRock( U,whichWay )
% convert veg & rock fractions within their total
% input
%   U - endmember structure with vegetation and rock fractions
%   whichWay - either 'forward' or 'inverse'
%       if 'forward', convert to fractions of their total
%       if 'inverse', convert back to fractions of total
% output
%   convS - same as input structure, except
%       if 'forward', rock and veg fractions sum to 1
%       if 'inverse', rock and veg fractions are their fraction of total

fillmember = {'vegetation_fraction','rock_fraction'};
whichWay = lower(whichWay);
convS = U;

% only do it if both rock and veg are present
if isfield(U,fillmember{1}) && isfield(U,fillmember{2})
    % make sure all values inside the mask are not NaN
    endmember = {'vegetation_fraction','rock_fraction',...
        'snow_fraction','other_fraction'};
    for k=1:length(endmember)
        if isfield(U,endmember{k})
            X = U.(endmember{k});
            X(isnan(X)) = 0;
            U.(endmember{k}) = setMask(X,U.imgMask);
        end
    end
    
    % forward or inverse
    switch whichWay
        case 'forward'
            % normalize the veg and rock fractions to sum to 1.0
            if isfield(U,'other_fraction')
                residual = 1-(U.snow_fraction+U.other_fraction);
            else
                residual = 1-U.snow_fraction;
            end
            t = residual>U.threshold;
            for m=1:length(fillmember)
                convS.(fillmember{m})(t) = U.(fillmember{m})(t)./residual(t);
                convS.(fillmember{m})(~t) = 0;
            end
            
        case 'inverse'
            % renormalize the vegetation and rock to their fraction of total
            sumT = U.(fillmember{1})+U.(fillmember{2}); % should = 1.0, except where all snow
            if isfield(U,'other_fraction')
                residual = 1-(U.snow_fraction+U.other_fraction);
            else
                residual = 1-U.snow_fraction;
            end
            residual(residual<0) = 0; % if exists, just rounding error
            t = residual>U.threshold & sumT>U.threshold;
            for m=1:length(fillmember)
                convS.(fillmember{m})(t) = (U.(fillmember{m})(t)./sumT(t)).*residual(t);
                convS.(fillmember{m})(~t) = 0;
            end
            
        otherwise
            error('argument %s not recognized',whichWay)
    end
end
end