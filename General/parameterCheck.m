function parameterCheck(Nrequired,functionName,varargin)
%check to see that all optional arguments are parameters
%
if isempty(varargin)
    return
end
Nopt = length(varargin);
if mod(Nopt,2)~=0
    warning('number of optional arguments must be even, checking for the error')
end
for k=1:2:Nopt
    assert(ischar(varargin{k}) || isstring(varargin{k}),...
        'in function %s, all optional arguments are name-value pairs, argument #%d is of class %s',...
        functionName,Nrequired+k,class(varargin{k}))
    assert(k+1<=Nopt,'in function %s, argument %d (''%s'') must be followed by its value',...
        functionName,Nrequired+k,varargin{k})
end
return
end