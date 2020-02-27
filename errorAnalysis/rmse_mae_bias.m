function [RMSE,MAE,bias] = rmse_mae_bias(yHat,y)
%RMSE Summary of this function goes here
%   yHat - modeled variable, algorithm output
%   y - truth measurement/observation

validateattributes(yHat,{'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', ...
    'single','double'},{'nonsparse'},mfilename,'A',1);
validateattributes(y,{'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', ...
    'single','double'},{'nonsparse'},mfilename,'B',1);

if ~isa(yHat,class(y))
    error(message('images:validate:differentClassMatrices','A','B'));
end
    
if ~isequal(size(yHat),size(y))
    error(message('images:validate:unequalSizeMatrices','A','B'));
end

if isempty(yHat) % If x is empty, y must also be empty
    RMSE = [];
    MAE = [];
    bias = [];
    return;
end

if isinteger(yHat)     
    yHat = double(yHat);
    y = double(y);
end

sumDiff = norm(yHat(:)-y(:),2);
absDiff = norm(yHat(:)-y(:),1);
N = numel(yHat);
MSE = sumDiff.^2/N;
RMSE = sqrt(MSE);
MAE= absDiff/N;
bias = sumDiff/N;
end

