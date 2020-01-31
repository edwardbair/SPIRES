function [ filteredCube ] = filterCloudCube( cloudNotSnow,maxDays,maxpctile )
% [ filteredCube ] = filterCloudCube( cloudNotSnow,maxDays,maxpctile )
%Filter multi-day cloud cube to eliminate clouds that are too persistent,
%specified both as maxDays (assumed to be snow instead) and maxpctile
%(upper percentile of consecutive days, also assumed to be snow).
%Threshold value used is the mean of the two.
%
% Input
%   cloudNotSnow - output from MODISclouds
%   maxDays - maximum consecutive days (set to zero to use maxpctile only)
%   maxpctile - upper percentile of the distribution, above which clouds
%       are eliminated (set to zero to use maxDays only)
%
% Output
%   filteredCube - input cloud cube filtered to eliminate spurious clouds

% median filter first, but only subtract clouds, don't add
C = false(size(cloudNotSnow));
parfor k=1:size(C,3)
    C(:,:,k) = medfilt2(cloudNotSnow(:,:,k),'symmetric');
end
t = C>cloudNotSnow;
C(t) = false;

% return if number of days in cube is less than maxDays
if size(C,3)<maxDays
    filteredCube = C;
    return
end

% image of total cloudy days for each pixel
C = reshape(C,size(C,1)*size(C,2),size(C,3))';
X = sum(C,1);
% maximum consecutive cloudy
maxConsecutive = ones(size(X)); % ok to lump zero and one
parfor k=1:size(C,2)
    if X(k)>1
        m = true(1,X(k)); % maximum possible consecutive for this pixel;
        v = C(:,k)'; % make this pixel a row vector
        mlen = length(m);
        while isempty(strfind(v,m(1:mlen)))
            mlen = mlen-1;
        end
        maxConsecutive(k) = mlen;
    end
end

% define threshold, based on percentiles and max expected
assert(maxDays>0 || maxpctile>0,...
    'must specify either maxDays or maxpctile, or both')
if maxDays==0
    thresh = prctile(maxConsecutive,maxpctile);
elseif maxpctile==0
    thresh = maxDays;
else %usual case
    thresh = mean([maxDays prctile(maxConsecutive,maxpctile)]);
end
% find days that exceed that threshold
t = maxConsecutive>=thresh;
m = true(1,ceil(thresh));
parfor k=1:length(t)
    if t(k)
        v = C(:,k)'; % make this pixel a row vector
        p = strfind(v,m);
        while ~isempty(p)
            % shorten the cloudy day period
            v(p(1)+length(m)-1) = false;
            p = strfind(v,m);
        end
        C(:,k) = v'; % back to column vector and replace in C
    end
end

filteredCube = reshape(C',size(cloudNotSnow));

end