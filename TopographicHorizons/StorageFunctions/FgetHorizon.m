function [horizon] = FgetHorizon(horizonInterpolant,row,col,azimuth)
%get the horizon for a set of rows, columns, and azimuths
%
%Input
% horizonFunction - output of SaveHorizon with the .mat option
% row - rows of the horizon of dimension (nRows,nCols,nAzimuths)
% col - columns of the horizon
% azimuth - azimuths of the horizon
% (row, col, and azimuth must be of the same dimension if not scalars)
%Output
% horizon - horizon angles (elevation from horizontal) for the inputs
%
%Note
% azimuths created by SaveHorizon and the inputs to this function must
% correspond to the convention specified in azimuthPreference (either -/+
% 180 degrees or 0 to 360)

%check sizes
[r,c,a] = checkSizes(row,col,azimuth);
%output
horizon = horizonInterpolant(r,c,a);
end