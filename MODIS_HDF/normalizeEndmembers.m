function N = normalizeEndmembers(F, endmember, varargin )
%Adjust endmember fractions to sum to 1.0
% N = normalizeEndmembers(F, endmember [, threshold] )
%
% Input
%   F - structure, with endmembers
%   endmember - cell vector of endmembers to include in the normalization
% Optional input
%   threshold - minimum value below which endmember fraction is set to zero
%
% Output
%   N- structure with normalized endmember planes or cubes, in same order
%       as endmember vector
%       (endmembers not normalized are not included)

% set threshold
numArg = 2;
if nargin>numArg
    threshold = varargin{1};
else
    threshold = 0;
end

% endmember names and sizes
siz = size(F.(endmember{1}));
numEnd = length(endmember);
for k=2:numEnd
    assert(isequal(size(F.(endmember{k})),siz),...
        'size of endmember %s different than endmember %s',...
        endmember{k},endmember{1})
end

% no endmember can exceed 1.0 or be negative
for k=1:numEnd
    t = F.(endmember{k})>1;
    if nnz(t)
        F.(endmember{k})(t) = 1;
    end
    t = F.(endmember{k})<threshold;
    if nnz(t)
        F.(endmember{k})(t) = 0;
    end
end

% output endmembers
for k=1:numEnd
    N.(endmember{k}) = zeros(siz,'like',F.(endmember{k}));
end

% normalize, each time plane independently
for s=1:siz(3)
    sumEnd = zeros(siz(1),siz(2));
    for k=1:numEnd
        tempN.(endmember{k}) = F.(endmember{k})(:,:,s);
        sumEnd = nansum(cat(3,sumEnd,tempN.(endmember{k})),3);
    end
    toAdjust = ~(sumEnd==1 | sumEnd==0);
    for k=1:numEnd
        tempN.(endmember{k})(toAdjust) =...
            tempN.(endmember{k})(toAdjust)./sumEnd(toAdjust);
        N.(endmember{k})(:,:,s) = tempN.(endmember{k});
    end
end
end