function [ output,varargout ] = interpAvg1( x,y,varargin )
% yq = interpAvg1(x,y,xq, [options] )
% pp = interpAvg1(x,y,'pp', [options])
%Interpolate averaged data (i.e. step function) where the abscissa values
%are at the end of each averaging period.
%This is often the case with meteorological data, where the instruments
%average over some period (e.g., hour) and report the average at the end of
%the time step.
%The same approach, with modification, can also be used to estimate an
%empirical PDF.
%
%Input
%   x - values of the abscissa, must continuously increase but need not be
%       equally spaced
%   y - values on the ordinate, y(n) is the average from x(n-1) to x(n)
%   xq - abscissa values at which the interpolated values are wanted
%   'pp' - result is returned as spline coefficients in pp format
%
%Optional input, name-value pairs
%   'firstDeltaX' - x(1)-x(0): y(1) represents the average over the time
%       step prior to the beginning of the sequence (default is same as
%       x(2)-x(1))
%   'smooth' - then next argument is either
%       'positive' - all values are non-negative
%       'negative' - all values are non-positive
%       'either' - values can be positive or negative

p = inputParser;
defaultXQ = [];
defaultDeltaX = [];
validationFcn = @(x) isnumeric(x) || ischar(x);
addRequired(p,'x',@isnumeric);
addRequired(p,'y',@isnumeric);
addOptional(p,'xq',defaultXQ,validationFcn);
addParameter(p,'deltaX',defaultDeltaX,@isnumeric);
addParameter(p,'smooth','',@ischar);
parse(p,x,y,varargin{:});

%% check sizes and constraints
assert(isrow(x) || iscolumn(x),'x and y must vectors')
rowform = isrow(x); % check if rows, reset at end
x = x(:);
y = y(:);
assert(length(x)>=3,'length(x) must be >=3');
assert(isequal(size(x),size(y)),'x and y must have same size');

% get rid of NaNs or infinite values
if any(isnan(y)|isinf(y))
    x(isnan(y)|isinf(y)) = [];
    y(isnan(y)|isinf(y)) = [];
    warning('NaNs or +/-Inf in input y ignored')
    assert(length(x)>=3,'length(x) must be >=3 after removal of NaNs');
end

% strip out duplicates
if ~isequal(x,unique(x))
    [C,ia,~] = unique(x);
    kd = find(diff(ia)~=1);
    newy = y(ia);
    for k=1:length(kd)
        kdup = find(x==C(kd(k)));
        if ~isempty(kdup)
            newy(kd(k)) = mean(y(kdup));
        end
    end
    x = C;
    y = newy;
    warning('duplicate x-values found, averaging y-values')
end

%% actual values or pp form
if isnumeric(p.Results.xq)
    xq = p.Results.xq(:);
    ppform = false;
    if ~isequal(xq,unique(xq))
        warning('xq is being resorted into increasing order, with no duplicates')
        xq(isnan(xq)|isinf(xq)) = [];
        xq = unique(xq);
    end
    if max(xq)>max(x)
        warning('xq truncated beyond last legitimate x value')
        xq(xq>max(x)) = [];
    end
    if nargout>1
        varargout{1} = xq;
    end
else % must be character string
    assert(~isempty(p.Results.xq),'3rd argument must be a vector or ''pp''');
    ppform = true;
    assert(strcmpi(p.Results.xq,'pp'),...
        'if 3rd argument not a vector of abscissa values, must be ''pp''')
end

%% smooth or not
if isempty(p.Results.smooth)
    doSmooth = false;
else
    doSmooth = true;
    prescription = slmset('interiorknots','fixed','knots',ceil(length(x)*2/3));
    switch p.Results.smooth
        case 'positive'
            prescription = slmset(prescription,'increasing','on','leftminvalue',0);
        case 'negative'
            prescription = slmset(prescription,'decreasing','on','leftmaxvalue',0);
        case 'either'
        otherwise
            error('''smooth'' value (''%s'') not recognized',p.Results.smooth)
    end
end

%% fit cumulative sum, then differentiate
% make dx same length as x
dx = diff(x);
if isempty(p.Results.deltaX)
    dx = cat(1,x(2)-x(1),dx);
else
    assert(p.Results.deltaX>0,'firstDeltaX, if specified, must be positive')
    dx = cat(1,deltaX,dx);
end

if doSmooth
    if ppform
        prescription = slmset(prescription,'result','pp');
        pp = slmengine(x,cumsum(y.*dx),prescription);
        output = fnder(pp);
    else
        slm = slmengine(x,cumsum(y.*dx),prescription);
        output = slmeval(xq,slm,1);
        if rowform && iscolumn(output)
            output = output';
        end
    end
    
else % not smoothed, just use pchip (this is best when data are accurate)
    Fp = fit(x,cumsum(y.*dx),'pchip');
    if ppform
        output = fnder(Fp.p);
    else
        output = differentiate(Fp,xq);
        if rowform && iscolumn(output)
            output = output';
        end
    end
end
end