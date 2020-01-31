function [ X, varargout ] = GetTopography( h5file, whichVariable )
% X = GetTopography( h5file, whichVariable )
% [ X, hdr ] = GetTopography( h5file, whichVariable )
%
% retrieve elevation, slope, aspect, or view factor from HDF 5 topography file
%
% Input
%   h5file - HDF 5 filename
%
% Output
%   X - output matrix of requested variable, converted to single-precision
%       floating point
% Optional output
%   hdr - structure with coordinate information about the dataset

% check to make sure input is a .h5 file
assert(strcmpi('.h5',h5file(end-2:end)),'input must be HDF 5 file')

% check if variable supported, case insensitive, only first 4 characters
Var = lower(whichVariable(1:4));
assert(strcmp(Var,'elev') ||...
    strcmp(Var,'slop') ||...
    strcmp(Var,'aspe') ||...
    strcmp(Var,'view'),...
    'input variable ''%s'' not recognized',whichVariable)
switch Var
    case 'elev'
        whichVariable = 'elevation';
    case 'slop'
        whichVariable = 'slope';
    case 'aspe'
        whichVariable = 'aspect';
    case 'view'
        whichVariable = 'viewfactor';
end
dataset = ['/Grid/' whichVariable];

% input data, converted to single if necessary
data = h5read(h5file,dataset);
I = h5info(h5file,dataset);
t = false(size(data));
if I.FillValue~=0
    t = data==I.FillValue;
else
    is16bitSignedInteger = strcmpi(I.Datatype.Class,'H5T_INTEGER') &&...
        strcmpi(I.Datatype.Type,'H5T_STD_I16LE') &&...
        I.Datatype.Size==2;
    if is16bitSignedInteger
        t = data==intmin('int16');
    end
end
divisor = h5readatt(h5file,dataset,'divisor');
X = single(data)/divisor;
if nnz(t)
    X(t) = NaN;
end

% other output argument?
nout = max(nargout,1) - 1;
if nout==1
    hdr = GetCoordinateInfo(h5file,'/Grid',size(X));
    varargout{1} = hdr;
end
end