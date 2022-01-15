function [V,method,catValues,catIndex,fillV] = fwdCategorical(X,method,fillvalue)
% convert categorical variable to signed integer but keep track of
% categories

assert(iscategorical(X),'input variable must be categorical')
if ~strcmpi(method,'nearest')
    warning('for categorical variable, interpolation method changed to ''nearest''')
    method = 'nearest';
end

%unique categories
catValues = categorical(categories(X(:)));
if ~isempty(fillvalue) && ~iscategorical(fillvalue)
    error('if ''fillvalue'' is specified for categorical input, it must also be categorical')
end
if length(catValues)<=intmax('int8') % probable answer
    catType = 'int8';
elseif length(catValues)<=intmax('int16') % unlikely
    catType = 'int16';
else
    error('input has %d distinct categories? really unlikely',...
        length(catValues))
end

catIndex = cast(1:length(catValues),catType);
% assign fillvalue index
if ~isempty(fillvalue)
    % is fillvalue among the categories
    t = catValues==fillvalue;
    if nnz(t)
        catIndex(t) = intmin(catType);
    else
        catIndex(end+1) = intmin(catType);
        catValues(end+1) = fillvalue;
    end
else
    catIndex(end+1) = intmin(catType);
    catValues(end+1) = categorical({'fill'});
end

V = zeros(size(X),catType)+intmin(catType);
for k=1:length(catValues)
    t = X==catValues(k);
    V(t) = catIndex(k);
end

fillV = intmin(catType);

end