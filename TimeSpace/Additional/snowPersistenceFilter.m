function Z = snowPersistenceFilter(Y,nPersist,threshold)
% Z = snowPersistenceFilter(Y,nPersist,threshold)
% remove days that persist for sequences shorter than nPersist
N = size(Y);
% reshape Y to put sequence of days in proximity in
% memory, as a matrix of size [N(3) N(1)*N(2)]
Y = reshape(Y,N(1)*N(2),N(3))';

% count NaNs as zero in the analysis, but restore them
tnan = isnan(Y);
Y(tnan) = 0;
% set to zero all columns whos total days are less than nPersist
isEndmember = Y>=threshold;
tot = sum(isEndmember,1);
t = tot<nPersist;
Y(:,t) = 0;

% only consider the remaining columns
t = ~t;
% for col=find(t), but using parfor
parfor col=1:length(t)
    if t(col)       
        [~,thisColumn] = persist1D(Y(:,col),nPersist,threshold,'minimum');
        Y(:,col) = thisColumn;
    end
end
% restore NaNs and reshape result back to original cube
if ~islogical(Y)
Y(tnan) = NaN;
end
Z = reshape(Y',N(1),N(2),N(3));
end