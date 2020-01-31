function U = float2integer(X,divisor,offset,inttype,varargin)
% U = float2integer(X,divisor,offset,inttype,varargin)
% convert float to scaled integer of specified type (to make files smaller)
% such that X can be recovered by X = offset + U/divisor
%
% input
%   X - floating point value (single or double) to be converted
%       by X = offset + U/divisor
%   divisor
%`  offset (in most cases, offset=0)
%   inttype - e.g. 'uint16', 'int16', 'uint8', or 'int8'
%
% optional input
%   minimum and maximum - either as 2 variables, or as vector of length 2
%
% Example -- scale aspects, whos values range from -180 to +180, as int8
%   I = float2integer(aspect,intmax('int8')/180,0,'int8',-180,180)

noptarg = size(varargin,2);
if noptarg>0
    if noptarg==1
        minmaxV = varargin{1};
        assert(length(minmaxV)==2,'if 1 optional arg, must be vector of length 2')
    elseif noptarg==2
        minmaxV = [varargin{1} varargin{2}];
        assert(length(minmaxV)==2,'if 2 optional args, must be scalars')
    else
        error('just 1 or 2 optional arguments, no more!')
    end
    X = truncateLimits(X,minmaxV);
end

X = double(X);
t = isnan(X);
% make sure stays within bounds of this integer type
limits = offset+[double(intmin(inttype))/divisor double(intmax(inttype))/divisor];
tt = ~t(:);
if any(tt)
    assert(nanmin(X(:))>=limits(1) && nanmax(X(:))<=limits(2),...
        'after conversion, some X values (%f<%f or %f>%f) fall outside integer %s''s range [%d %d]',...
        nanmin(X(:)),limits(1),nanmax(X(:)),limits(2),intmin(inttype),intmax(inttype))
end

U = cast(round((X-offset)*divisor),inttype);
U(t) = intmin(inttype);

end