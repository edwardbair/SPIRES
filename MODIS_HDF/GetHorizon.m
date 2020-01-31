function [ X, varargout ] = GetHorizon( h5file, azm, varargin )
% angles = GetHorizon(h5file,azimuth(s))
% [angles,hdr] = GetHorizon(h5file,azimuth(s))
% mask = GetHorizon(h5file,azimuth(s),solarzenith(s))
% [mask,hdr] = GetHorizon(h5file,azimuth(s),solarzenith(s))
%
% retrieve horizon angles or horizon mask
%
% Input
%   h5file - HDF 5 filename
%   azimuth - in degrees, scalar, or matrix of size of topographic grid,
%       for the direction(s) in which the horizon angles are desired
%       [azimuth to south is 0 degrees, + east of south, - west of south]
% Optional input
%   zenith - solar zenith angle(s), scalar or matrix of size of topographic
%       grid, which if specified causes the output to be a mask, 1 for
%       cells that the sun illuminates, 0 for points that are hidden
%
% Output
%   if just azimuth input argument, angles to the horiz
%   if azimuth and zenith are given, horizon mask with illuminated points =
%       1 and hidden points 0
% Optional output
%   hdr - structure with coordinate information about the dataset

% check to make sure input is a .h5 file
assert(strcmpi('.h5',h5file(end-2:end)),'input must be HDF 5 file')
info = h5info(h5file);
siz = info.Groups(1).Datasets(1).ChunkSize;
if ndims(siz)==2
    siz = [siz(1) siz(2) 1];
end

% heading information
azm = round(azm); % just need to the nearest degree
u = unique(azm(~isnan(azm)));
nHorz = h5readatt(h5file,'/Grid/horizons','nHorizons');
hbound = h5readatt(h5file,'/Grid/horizons','AzmLimits');
divisor = h5readatt(h5file,'/Grid/horizons','divisor');
angles = linspace(hbound(1),hbound(2),nHorz);
zeroindex = find(angles==0);
uindex = u+zeroindex;

% solar zenith too
nVarargs = length(varargin);
doMask = false;
if nVarargs==1
    solarZ = varargin{1};
    doMask = true;
    cosZ = cosd(solarZ);
end

% sines of the horizon angles
if isempty(u) % entire grid in darkness, so just need to get size
    iH = h5read(h5file,'/Grid/horizons',[1 1 1],siz);
    if doMask
        X = false(size(iH));
    else
        X = ones(size(iH),'single'); % set horizon to 90 deg where sun not up
    end
elseif length(u)==1
    iH = h5read(h5file,'/Grid/horizons',[1 1 uindex],siz);
    sH = double(iH)/divisor;
else
    for k=1:length(uindex)
        iH = h5read(h5file,'/Grid/horizons', [1 1 uindex(k)], siz);
        if k==1
            output = zeros(size(iH),'like',iH);
            assert(isscalar(azm) || isequal(size(azm),size(iH)),...
                'if azimuth is not a scalar, it must be size of horizon grid')
        end
        t = azm==u(k);
        output(t) = iH(t);
    end
    % all of output should be filled, except for cells where sun has not
    % yet risen
    sH = double(output)/divisor;
    sH(isnan(output)) = 1; % set horizon to 90 deg for cells where sun has not risen
end

% mask, or angles themselves?
if ~isempty(u)
    if doMask
        assert(isscalar(cosZ) || isequal(size(cosZ),size(iH)),...
            'if solarZ is not a scalar, it must be size of horizon grid')
        X = sH<cosZ;
        X(cosZ<=0 | isnan(cosZ) | isnan(azm)) = false;
    else
        X = single(asind(sH));
    end
end

% other output argument?
nout = max(nargout,1) - 1;
if nout==1
    hdr = GetCoordinateInfo(h5file,'/Grid',size(X));
    varargout{1} = hdr;
end
end