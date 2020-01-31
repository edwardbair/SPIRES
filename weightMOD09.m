function [ matdate,weights,bandWeights ] = weightMOD09( MOD09filelist,varargin )
% [ matdate,weights,bandWeights ] = weightMOD09( MOD09filelist,[topographyFile] )
%examine multi-day set of MOD09 or MYD09 images (typically 8-32 days) and
%identify which pixels we have most confidence in on a scale from 0 to 1
%
% (some analysis based on Karl Rittger's modis8dayPixelRank.m
%
% Input
%   filelist - list of MODIS MOD09 of MYD09 files, generally between 8 and 32
%       days (to cover 1/2, 1, or 2 16-day repeat cycles), but can be just
%       1 day or can be whole year
% Optional input
%   topographyFile - our H5 topography file corresponding to the MODIS tile
%
% Output
% vector size [nDay 1]
%   matdates - MATLAB datenums of the files
%   weights - array of size (M x N x nDay), confidence weights (0-1) based
%       on likely quality of each pixel, which can be used for selecting,
%       smoothing, or filtering (bad pixels are set to weight=0)
%   bandWeights - cell variable of length nDay, with each cell an array of
%       size (M x N x nBands) with weights for each band
%
% Jeff Dozier, 2014-05-25

% dates of the files, also make sure the prefix and tile are the same
if iscell(MOD09filelist)
    nFiles = length(MOD09filelist);
else
    nFiles = 1;
    MOD09filelist = {MOD09filelist};
end
fT = infoFromMODISfilename(MOD09filelist);
matdate = fT{:,'matdate'};
for k=2:nFiles
    assert(strcmpi(fT{k,'prefix'},fT{1,'prefix'}),...
        'prefix %d (%s) different than prefix 1 (%s)',...
        k,fT{k,'prefix'},fT{1,'prefix'})
    assert(strcmpi(fT{k,'tile'},fT{1,'tile'}),...
        'tile %d (%s) different than tile 1 (%s)',...
        k,fT{k,'tile'},fT{1,'tile'})
end

% slope & aspect from topography file, if present
if nargin>1
    useTopo = true;
    slope = GetTopography(varargin{1},'slope');
    aspect = GetTopography(varargin{1},'aspect');
else
    useTopo = false;
end

for d=1:height(fT)
    weight1km = GetWeight1km(fT{d,'file'});
    if useTopo
        weightAngles = GetWeightAngles(fT{d,'file'},slope,aspect,varargin{1});
    else
        weightAngles = GetWeightAngles(fT{d,'file'});
    end
    weight500m = GetWeight500m(fT{d,'file'});
    if useTopo
        W = cat(3,imresize(weight1km,2),weightAngles,weight500m);
    else
        W = cat(3,imresize(weight1km,2),imresize(weightAngles,2),weight500m);
    end
    weightBands = GetWeightBands(fT{d,'file'});
    if d==1
        weights = nanmin(W,[],3);
        bandWeights = cell(size(fT,1),1);
    else
        weights = cat(3,weights,nanmin(W,[],3));
    end
    weights(weights<0) = 0;
    weightBands(weightBands<0) = 0;
    bandWeights{d} = weightBands;
end

end

function w = assignWeights(stateVariable,values,weights)
w = ones(size(stateVariable),'single');
for k=1:length(values)
    t = stateVariable==values(k);
    w(t) = weights(k);
end
end

function W = GetWeight1km(file)
% assign weights for the 1 km state data
S = GetMOD09GA(file,'state');
siz = size(S.cloud);

% create cube of weights corresponding to possible values of the state
% variable
% cloud flag
values = cast(0:3,'like',S.cloud);
weights = [1 .5 .5 1]; % assign 0.5 to cloudy, because we're not sure
W = assignWeights(S.cloud,values,weights);

% cloud shadow
shadowWeight = .8;
thisWeight = ones(siz,'single');
thisWeight(S.cloudshadow) = shadowWeight;
W = cat(3,W,thisWeight);

% land water
values = cast(0:7,'like',S.landwater);
weights = [.8 1 .8 .8 .8 0 0 0];
thisWeight = assignWeights(S.landwater,values,weights);
W = cat(3,W,thisWeight);

% aerosol
values = cast(0:3,'like',S.aerosol);
weights = [1 1 .5 .25];
thisWeight = assignWeights(S.aerosol,values,weights);
W = cat(3,W,thisWeight);

% cirrus
values = cast(0:3,'like',S.cirrus);
weights = [1 1 .5 .25];
thisWeight = assignWeights(S.cirrus,values,weights);
W = cat(3,W,thisWeight);

% snow - ignored

end

function w = GetWeightAngles(file,varargin)
earthRadius = 6.371007181e+03;
orbitHeight = 705;

useTopo = nargin>1;
if useTopo
    slope = varargin{1};
    aspect = varargin{2};
    tfile = varargin{3};
end

% sensor zenith, and pixel sizes
sensorZ = GetMOD09GA(file,'sensorzenith');
if useTopo
    sensorZ = imresize(sensorZ,2);
end
[ppl,ppt,~] = pixelSize(earthRadius,orbitHeight,1,sensorZ);
w = 1./(ppl.*ppt);

% solar zenith
solarZ = double(GetMOD09GA(file,'solarzenith'));
mu0 = cosd(solarZ);
if useTopo
    mu0 = imresize(mu0,2);
end

% topographic
if useTopo
    solarAz = 180-double(GetMOD09GA(file,'solarazimuth')); % convert to 0 south
    X = imresize(cat(3,sind(solarAz),cosd(solarAz)),2);
    solarAz = atan2d(X(:,:,1),X(:,:,2)); % atan2d(sinA,cosA)
    mu = sunslope( mu0,solarAz,slope,aspect );
    solarZ = imresize(solarZ,2);
    visible = GetHorizon(tfile,solarAz,solarZ);
    mu(~visible) = 0;
    w = cat(3,w,sqrt(mu));
else
    w = cat(3,w,sqrt(mu0));
end
end

function w = GetWeight500m(file)
% assign weights based on the 500m quality information

% modland quality flags
QC_500m = GetMOD09GA(file,'QC');
values = cast(0:3,'like',QC_500m.modland);
weights = [1 .5 0 0];
w = assignWeights(QC_500m.modland,values,weights);

% atmospheric correction
thisWeight = ones(size(QC_500m.atmosCorr),'single');
thisWeight(~QC_500m.atmosCorr) = 0.8;
w = cat(3,w,thisWeight);

% adjacency correction
thisWeight = ones(size(QC_500m.adjCorr),'single');
thisWeight(~QC_500m.adjCorr) = 0.9;
w = cat(3,w,thisWeight);

end

function w = GetWeightBands(file)
% weights assigned to each band
QC_500m = GetMOD09GA(file,'QC');
values = cast([0 7:15],'like',QC_500m.bandQA);
weights = [1 0 0 cosd(87) cosd(85.5) 0 .5 .25 0 0];
w = assignWeights(QC_500m.bandQA,values,weights);
end