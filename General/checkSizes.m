function [ varargout ] = checkSizes( varargin )
% [ varargout ] = checkSizes( varargin )
%general function to check sizes of inputs
%all inputs that are not scalars must be same size
%all inputs that are scalars are expanded to same size as non-scalar input(s)

% check number of input and output arguments
narginchk(1,Inf)
nargoutchk(nargin,nargin)

% which ones not scalar (set to scalar if just replicated scalar)
scalar = false(size(varargin));
for k=1:length(varargin)
    scalar(k) = isscalar(varargin{k});
    if ~scalar(k)
        if isscalar(unique(varargin{k}(:)))
            scalar(k) = true;
            varargin{k} = unique(varargin{k}(:));
        end
    end
end

% trivial case, all scalar
if all(scalar)
    for k=1:length(varargin)
        varargout{k} = varargin{k}; %#ok<*AGROW>
    end
    
    % just one non-scalar, set others to same size
elseif sum(~scalar)==1
    n = find(~scalar);
    for k=1:length(varargin)
        if k==n
            varargout{k} = varargin{k};
        else
            varargout{k} = repmat(varargin{k},size(varargin{n}));
        end
    end
    
    % general case
else
    % vectors must be row or column, but not both
    for k=1:length(varargin)
        if isrow(varargin{k})
            varargin{k} = varargin{k}.'; % must use the .' as variable might be complex
        end
    end
    n = find(~scalar);
    N = size(varargin{n(1)});
    for k=2:length(n)
        assert(isequal(N,size(varargin{n(k)})),...
            'variables that are not scalars must have same size, variable %d has different size than variable %d',n(k),n(1))
    end
    for k=1:length(varargin)
        if scalar(k)
            varargout{k} = repmat(varargin{k},N);
        else
            varargout{k} = varargin{k};
        end
    end
end
end