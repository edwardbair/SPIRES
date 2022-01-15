function [id,varargout]=persist1D(x,nPersist,threshold,option)
% [id [,newx]]=persist1D(x,nPersist,minmaxThresh,option)
% identify sequences of variable x that meet threshold values, depending on
% the value of option, for either at least nPersist or no greater than
% nPersist
%
% Input
%   x - vector of values, either numeric or logical
%   nPersist - duration of values to keep (can be min or max)
%   threshold - values below this to consider zero (ignored if x is logical)
%   option - either 'minimum' (shorter durations are considered false) or
%       'maximum' (longer durations are considered false)
% Output
%   id - logical vector if values that meet the criteria
% Optional output
%   newx - = original x with x(~id) set to zero

% identify the start [false true] and stop [true false] of all potential sequences
% use strfind, which requires a row vector instead of a column vector
if iscolumn(x)
    inputColumns = true;
    x = x';
else
    inputColumns = false;
end

% potential values are any above threshold
if islogical(x)
    pot = x;
else
    pot = x>=threshold;
end
beg = strfind(pot,[false true]);
endd = strfind(pot,[true false]);
% for beg, we want the beginning of the sequence of 1s so get past the 0s
beg = 1+beg;
% and the whole series could start with 1 or end with 1
if pot(1)
    beg = cat(2,1,beg);
end
if pot(end)
    endd = cat(2,endd,length(pot));
end
% must have same number of beginnings and endings
assert(isequal(size(beg),size(endd)),...
    'something wrong, not same number of beginnings and endings')

switch option
    case 'minimum'
        % eliminate sequences that are not long enough
        longEnough = endd-beg>=nPersist-1;
        beg = beg(longEnough);
        endd = endd(longEnough);
    case 'maximum'
        % eliminate sequences that are too long
        shortEnough = endd-beg<=nPersist-1;
        beg = beg(shortEnough);
        endd = endd(shortEnough);
    otherwise
        error('option = %s not recognized',option)
end

% set the output logical variable to true between all beg-endd intervals
id = false(size(x));
for k=1:length(beg)
    id(beg(k):endd(k)) = true;
end

% switch back to column vector
if inputColumns
    id = id';
    x = x';
end
if nargout>1
    if logical(x)
        newx = false(size(x));
    else
        newx = zeros(size(x));
    end
    newx(id) = x(id);
    varargout{1} = newx;
end
end