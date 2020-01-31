function Z = cloudPersistenceFilter(Y,nPersist)
% remove cloudy days that persist for sequences longer than nPersist (in
% this case, Y is a logical variable)
assert(islogical(Y),'input Y must be logical')
N = size(Y);
% reshape Y to put sequence of days in proximity in
% memory, as a matrix of size [N(3) N(1)*N(2)]
Y = reshape(Y,N(1)*N(2),N(3))';

% only consider columns whose total cloudy days are greater than nPersist
tot = sum(Y,1);
t = tot>nPersist;

% for col=find(t), but using parfor
lengthColumn = size(Y,1);
parfor col=1:lengthColumn
    if t(col)
        id = persist1D(Y(:,col),nPersist,0,'maximum');
        Y(:,col) = id;
    end
end
% reshape result back to original cube
Z = reshape(Y',N(1),N(2),N(3));
end