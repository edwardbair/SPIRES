function y=TabIntegral(xt,yt,a,b)
% y=TabIntegral(xt,yt,a,b)
% compute integral from tabular data using integral + interp1
% input
%   xt,yt - input tabular data, sorted in xt
%   a,b - limits for integration

% "function" consists of values of yt at points xt (sorted)
y=integral(@tabfun,a,b);

    function s=tabfun(xi) % value at particular values of xi
        s=interp1(xt,yt,xi,'pchip');
    end
end