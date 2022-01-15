function [ X, varargout ] = h5getHorizon( h5file, varargin )
% angles = h5getHorizon(h5file)
% angles = h5getHorizon(h5file,azimuth)
% [___,hdr] = h5GetHorizon(h5file,____)
%
% retrieve horizon angles or horizon mask
%
% Input
%   h5file - HDF 5 filename
% Optional input, in order
%   azimuth - in degrees, either scalar or matrix of size of topographic grid,
%       for the direction(s) in which the horizon angles are desired
%       (orientation set by azimuthPreference, either azimuth to south is 0
%       degrees, + east of south, - west of south, or azimuth clockwise
%       from north
%
% Output
%   if just azimuth input argument, angles to the horizon in that direction
%   otherwise, horizon angles in all directions
% Optional output
%   hdr - structure with coordinate information about the dataset

p = inputParser;
addRequired(p,'h5file',@ischar)
addOptional(p,'azimuth',[],@isnumeric)
parse(p,h5file,varargin{:});
azm = p.Results.azimuth;

% check to make sure input is a .h5 file
[~,~,ext] = fileparts(h5file);
assert(strcmpi(ext,'.h5'),'input must be HDF 5 file')
info = h5info(h5file);
siz = info.Groups(1).Datasets(1).ChunkSize;
FillValue = info.Groups(1).Datasets(1).FillValue;
if ismatrix(siz)
    siz = [siz(1) siz(2) 1];
end

% heading information
angles = h5readatt(h5file,'/Grid','azimuths');
divisor = h5readatt(h5file,'/Grid/horizons','Divisor');

% horizon angles, get all if input azimuth empty
if isempty(azm)
    iH = h5read(h5file,'/Grid/horizons');
    X = double(iH)/divisor;
elseif isscalar(azm)
    adiff = abs(angles-azm);
    idx = find(adiff==min(adiff));
    iH = h5read(h5file,'/Grid/horizons',[1 1 idx],siz);
    X = double(iH)/divisor;
    X(iH==FillValue) = NaN;
else
    assert(isequal(size(azm),siz(1:2)),...
        'if azimuth is not a scalar, it must be size of horizon grid')
    u = unique(azm(:));
    for k=1:length(u)
        adiff = abs(angles-u(k));
        idx = find(adiff==min(adiff));
        iH = h5read(h5file,'/Grid/horizons', [1 1 idx], siz);
        if k==1
            output = zeros(size(iH),'like',iH);
        end
        t = azm==u(k);
        output(t) = iH(t);
    end
    X = double(output)/divisor;
    % all of output should be filled, except for cells where sun has not
    % yet risen
    X(output==FillValue) = NaN;
end

% other output argument?
if nargout>1
    hdr = h5getCoordinates(h5file);
    varargout{1} = hdr;
end
end