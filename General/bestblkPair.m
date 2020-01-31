function [ blksiz ] = bestblkPair( matrixSize, nProc, varargin )
%find pair of block sizes that have no overlapping seams
% [ blksiz ] = bestblkPair( matrixSize, nProc [, maxval] )
%   use in block processing where edge effects might exist
% Input
%   matrixSize - dimensions of matrix, vector of length 2
%   nProc - number of processors, can be either number of nodes or number of
%       cores, depending on how processing is distributed
% Optional input
%   maxval - maximum block value
%
% Output
%   blksiz - 2x2 matrix, with each row containing a different blocksize

assert(numel(matrixSize)==2,'matrixSize must be a vector of 2 elements')

N = sqrt(nProc);
tooDifferent = true;
if nargin>2
    maxval = varargin{1};
else
    maxval = Inf;
end
while tooDifferent
    if matrixSize(1)==matrixSize(2)
        b = best1Dp(matrixSize(1),N,maxval);
        blksiz = [b b];
    else
        b1 = best1Dp(matrixSize(1),N,maxval);
        b2 = best1Dp(matrixSize(2),N,maxval);
        blksiz = [b1 b2];
    end
    r = max(blksiz(1,:)./blksiz(2,:));
    tooDifferent = r>3;
    N = N+1;
end

end

function b=best1Dp(nVals,N,maxval)
if nVals/N>maxval
    N = nVals/maxval;
end
overlap = true;
M = N+1/2;
b1 = bestblk([nVals 1],nVals/N);
b1(2) = [];
b2 = b1;
while overlap && b2>=b1/2
    b2 = bestblk([nVals 1],nVals/M);
    b2(2) = [];
    % check overlap
    seq1 = 1:b1:nVals;
    seq2 = 1:b2:nVals;
    seq1(1) = [];
    seq2(1) = [];
    overlap = false;
    for k=1:length(seq1)
        kf = seq2==seq1(k);
        if nnz(kf)
            overlap = true;
            break;
        end
    end
    M = M+1/2;
end
if b2<b1/2 % didn't find within acceptable range, so drive with different nVals
    M = N+1/2;
    b2 = b1;
    while overlap && b2>=b1/2
        b2 = bestblk([nVals+1 1],nVals/M);
        b2(2) = [];
        % check overlap
        seq1 = 1:b1:nVals;
        seq2 = 1:b2:nVals;
        seq1(1) = [];
        seq2(1) = [];
        overlap = false;
        for k=1:length(seq1)
            kf = seq2==seq1(k);
            if nnz(kf)
                overlap = true;
                break;
            end
        end
        M = M+1/2;
    end
end
b = [b1 b2]';
end