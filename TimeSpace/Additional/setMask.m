function mv = setMask(v,m,varargin)
% mv = setMask(v,m,varargin)
% optional argument is fill value, default NaN
% v is input matrix or cube, m is mask matrix
p = inputParser;
validationFcn = @(x) isnumeric(x) || islogical(x);
addRequired(p,'v',validationFcn)
addRequired(p,'m',@islogical)
addOptional(p,'setValue',[],@isnumeric)
parse(p,v,m,varargin{:});
% make sure sizes compatible, mask could be a cube or a plane
sizeV = size(v);
sizeM = size(m);
% optional argument for fill value
if isempty(p.Results.setValue)
    if isinteger(v)
        setValue = intmin(class(v));
    elseif islogical(v)
        setValue = false;
    else % v is a float
        setValue = NaN;
    end
else
    setValue = p.Results.setValue;
end

if isequal(sizeV,sizeM)
    v(~m) = setValue;
    mv = v;
elseif isequal(sizeV(1:2),sizeM)
    m = reshape(m,sizeM(1)*sizeM(2),1);
    v = reshape(v,sizeV(1)*sizeV(2),sizeV(3));
    for k=1:size(v,2)
        v(~m,k) = setValue;
    end
    mv = reshape(v,sizeV);
else
    error(['sizes of image (' num2str(sizeV) ') and mask ('...
        num2str(sizeM) ') must be same'])
end
end