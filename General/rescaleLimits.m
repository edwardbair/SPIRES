function [ Y ] = rescaleLimits( X,limits )
%rescale X to stay between limits
%
% Input
%   X - input data, vector or matrix or N-D
%   limits - vector of length 2: min and max of output
%
% Output
%   Y - rescaled X with min and max of limits
%       (NaNs in input are NaN in output)

y0 = min(limits);
y1 = max(limits);
x0 = nanmin(X(:));
x1 = nanmax(X(:));

% y = a*x+b
a = (y1-y0)/(x1-x0);
b = (x1*y0 - x0*y1)/(x1-x0);
Y = a*X+b;

end