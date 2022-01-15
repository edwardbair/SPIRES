function snowR = prepCloudSnowFilter(sensor, varargin)
% snowR = prepCloudSnowFilter(sensor, ...)
% Preparation for running filterCloudSnow
%   Run this before calling filterCloudSnow, and pass the output from this
%   function to the snowR input table of that function.
%
% Input
%   sensor - 'modis' (same as 'modisterra'), 'modisaqua' 'landsat' (meaning L8), or 'viirs' case-insensitive
%
% Optional input - name, value pairs
%   'solarZ' - solar zenith angle(s) to run, default full range
%   'minradius' - minimum grain radius for snow in mm, scalar, default .04
%   'maxradius' - maximum grain radius for snow in mm, scalar, default 1.5
%   'deltavis' - factor to degrade reflectivity in minimum threshold visible
%       wavelengths because of dust, default 0.3
%   'dropletradius' - cloud droplet radius, mm, default 5e-3
%   'thickcloudwater' - thick cloud water equivalent, default 1 kg/m^2
%   'thincloudwater' - thin cloud water equivalent, default 0.1
%   'cirruswater' - cirrus cloud water equivalent, default 0.05
%
% Output
%   snowR - snow reflectivity structure, contains angles,
%       and the following reflectivities for Vis,
%           NIR (up to 1100 nm), and SWIR bands, along with their min (snow)
%           or max (cloud) NDSI values
%       maximum clean fine snow
%       minimum clear coarse snow
%       minimum coarse dirty snow
%           (note: not dealing with shallow snow)
%       thick clouds
%       thin clouds
%       cirrus clouds

cosineMultiplier = 100;

% parse input
p = inputParser;
defaultMinRadius = .04; % mm
defaultMaxRadius = 1.5; % mm
% use fplot to get the right spacing for angles
fp = fplot(@cosd,[0 89]);
defaultSolarZ = fp.XData;
defaultDeltaVis = 0.3;
% cloud size stuff from
% http://www-das.uwyo.edu/~geerts/cwx/notes/chap08/moist_cloud.html
% and
% http://neo.sci.gsfc.nasa.gov/view.php?datasetId=MYDAL2_D_CLD_WP
defaultCloudDropletRadius = 5e-3; % mm
defaultCirrusRadius = min(90e-3,defaultMinRadius*3/5);
defaultThickCloudWaterEquiv = 1; % kg/m^2 or mm
defaultThinCloudWaterEquiv = .1;
defaultCirrusWaterEquiv = .05;
% inputs
addRequired(p,'sensor',@ischar)
addParameter(p,'solarZ',defaultSolarZ,@isnumeric)
addParameter(p,'minradius',defaultMinRadius,@isnumeric)
addParameter(p,'maxradius',defaultMaxRadius,@isnumeric)
addParameter(p,'deltavis',defaultDeltaVis,@isnumeric)
addParameter(p,'dropletradius',defaultCloudDropletRadius,@isnumeric)
addParameter(p,'cirrusradius',defaultCirrusRadius,@isnumeric)
addParameter(p,'thickcloudwater',defaultThickCloudWaterEquiv,@isnumeric)
addParameter(p,'thincloudwater',defaultThinCloudWaterEquiv,@isnumeric)
addParameter(p,'cirruswater',defaultCirrusWaterEquiv,@isnumeric)
parse(p,sensor,varargin{:})

% variables
angles = p.Results.solarZ;
minRadius = p.Results.minradius;
maxRadius = p.Results.maxradius;
dirtfactor = p.Results.deltavis;
cloudDropletRadius = p.Results.dropletradius;
cirrusRadius = p.Results.cirrusradius;
thickCloudWaterEquiv = p.Results.thickcloudwater;
thinCloudWaterEquiv = p.Results.thincloudwater;
cirrusWaterEquiv = p.Results.cirruswater;

% spectral "land" bands (thermal not included)
switch lower(sensor)
    case {'modis','modisterra'}
        spectralBands = [620 670; 841 876; 459 479; 545 565; 1230 1250;...
            1628 1652; 2105 2155];
        visBands = [1 3 4];
        nirBands = 2;
        swirBands = [5 6 7];
        cirrusBand = [];
    case 'modisaqua'
        spectralBands = [620 670; 841 876; 459 479; 545 565; 1230 1250;...
            1628 1652; 2105 2155];
        visBands = [1 3 4];
        nirBands = 2;
        swirBands = [5 7];
        cirrusBand = [];
    case 'landsat'
        spectralBands = [433 453; 450 515; 525 600; 630 680; 845 885;...
            1560 1660; 2100 2300; 500 680; 1360 1390];
        visBands = [2 3 4 8];
        nirBands = 5;
        swirBands = [6 7];
        cirrusBand = 9;
    case 'viirs'
        % just the 'imagery' bands at 375 nadir resolution, not including
        % the mid-IR or thermal bands (adjust this later)
        spectralBands = [600 680; 846 886; 1580 1640];
        visBands = 1;
        nirBands = 2;
        swirBands = 3;
        cirrusBand = [];
    otherwise
        error('sensor %s not recognized',sensor)
end
doCirrus = ~isempty(cirrusBand);

if isrow(angles)
    angles = angles';
end
assert(~any(angles<0) && ~any(angles>=90),...
    'solar zenith angle(s) restricted to 0<=angle<90 degrees')
% round cosine index to nearest integer
cosIndex = (min(floor(cosineMultiplier*cosd(angles))):max(ceil(cosineMultiplier*cosd(angles))))';
cosAngles = cosIndex/cosineMultiplier;

% threshold reflectivity in each band for each angle
maximumSnow = reflectivity('snow',cosAngles,spectralBands,minRadius,Inf);
minimumSnow = reflectivity('snow',cosAngles,spectralBands,maxRadius,Inf);
% degrade minimum snow for dirt
dirtySnow = minimumSnow;
dirtySnow(:,visBands) = minimumSnow(:,visBands)*(1-dirtfactor);
dirtySnow(:,nirBands) = minimumSnow(:,nirBands)*(1-dirtfactor/3);

% clouds
thickCloud = reflectivity('cloud',cosAngles,spectralBands,...
    cloudDropletRadius,thickCloudWaterEquiv);
thinCloud = reflectivity('cloud',cosAngles,spectralBands,...
    cloudDropletRadius,thinCloudWaterEquiv);
cirrus = reflectivity('cirrus',cosAngles,spectralBands,...
    cirrusRadius,cirrusWaterEquiv);

% fill the snowR structure
snowR.sensor = sensor;
snowR.cosineMultiplier = cosineMultiplier;
snowR.spectralBands = spectralBands;
snowR.visBands = visBands;
snowR.nirBands = nirBands;
snowR.swirBands = swirBands;
snowR.cirrusBand = cirrusBand;
% table of means across visible, nir, and swir
snowR.remark = 'tables of means -> vis, nir, swir';
meanMaxSnow = [mean(maximumSnow(:,visBands),2) mean(maximumSnow(:,nirBands),2)...
    mean(maximumSnow(:,swirBands),2)];
meanMinSnow = [mean(minimumSnow(:,visBands),2) mean(minimumSnow(:,nirBands),2)...
    mean(minimumSnow(:,swirBands),2)];
meanThickCloud = [mean(thickCloud(:,visBands),2) mean(thickCloud(:,nirBands),2)...
    mean(thickCloud(:,swirBands),2)];
meanThinCloud = [mean(thinCloud(:,visBands),2) mean(thinCloud(:,nirBands),2)...
    mean(thinCloud(:,swirBands),2)];
meanDirtySnow = [mean(dirtySnow(:,visBands),2) mean(dirtySnow(:,nirBands),2)...
    mean(dirtySnow(:,swirBands),2)];
if doCirrus
    meanCirrus = [mean(cirrus(:,visBands),2) mean(cirrus(:,nirBands),2)...
        cirrus(:,cirrusBand)];
else
    meanCirrus = [mean(cirrus(:,visBands),2) mean(cirrus(:,nirBands),2)...
        mean(cirrus(:,swirBands),2)];
end

vnames = {'cosIndex','fineSnow','coarseSnow','dirtySnow','thickCloud',...
    'thinCloud','cirrus'};
snowR.meanTable = table(cosIndex,meanMaxSnow,meanMinSnow,meanDirtySnow,meanThickCloud,...
    meanThinCloud,meanCirrus,...
    'VariableNames',vnames);
snowR.wholeTable = table(cosIndex,maximumSnow,minimumSnow,dirtySnow,thickCloud,thinCloud,...
    cirrus,...
    'VariableNames',vnames);
% best NDSI values for discrimination
snowV = [2 3 4];
cloudV = [5 6];
cirrusV = 7;
[snowR.NDSI,snowR.ndsiBands] = findBestNDSI(snowR,vnames,snowV,cloudV);
if ~isempty(snowR.cirrusBand)
    [NDcirrus,cbands] = findBestNDSI(snowR,vnames,snowV,cirrusV);
    snowR.ndsiBands(3) = cbands(2);
    snowR.NDSI.cirrus = NDcirrus.cirrus;
end

deltavis = dirtfactor;
snowR.inputSpecs = table(minRadius,maxRadius,deltavis,...
    cloudDropletRadius,cirrusRadius,...
    thickCloudWaterEquiv,thinCloudWaterEquiv,cirrusWaterEquiv);

end

% reflectivity of snow or cloud in each band
function R=reflectivity(object,mu0,lambda,radius,waterEquivalent,varargin)
% object: 'snow', 'cloud' or 'cirrus'
% lambda: in nm
% radius: in mm
% waterEquivalent: in mm
% optional argument, reflectivity of surface under the cloud

switch lower(object)
    case 'snow'
        substance = 'ice';
    case 'cloud'
        substance = 'water';
    case 'cirrus'
        substance = 'ice';
    otherwise
        error('object must be ''snow'', ''cloud'' or ''cirrus''')
end

numArg = 5;
if nargin>numArg
    R0 = varargin{1};
else
    R0 = .3;
end

% Mie calculations
nBands = size(lambda,1);
% calculate at nWaves wavelengths for each band, then we can average
% the reflectivity across the bandpass
nWaves = 5;
wavesPerBand = zeros(nBands,nWaves);
for b=1:nBands
    wavesPerBand(b,:) = linspace(lambda(b,1),lambda(b,2),nWaves);
end
rIndex = RefractiveIndex(wavesPerBand,substance,'nm');
M = MieSphere(radius,'mm',wavesPerBand,'nm','refindex',rIndex);
% M = MieSphere(radius,'mm',wavesPerBand,'nm',rIndex);

% optical depth
if isinf(waterEquivalent)
    tau = inf(size(M.Qext));
else
    if strcmpi(object,'snow') || strcmpi(object,'cirrus')
        tau = tauSnow(convertLengthUnits(radius,'mm','m'),waterEquivalent,M.Qext);
    else
        tau = tauCloud(convertLengthUnits(radius,'mm','m'),waterEquivalent,M.Qext);
    end
end

% reflectivities for the wavelengths at each angle
R = zeros(length(mu0),nBands);
tmpRefl = zeros(size(M.omega,2),1);
omega = M.omega;
g = M.g;
for a=1:length(mu0)
    for k=1:size(omega,1)
        parfor m=1:size(omega,2)
            [tmpRefl(m),~,~,~] =...
                twostream(mu0(a),omega(k,m), g(k,m), tau(k,m),'R0',R0); %#ok<PFBNS>
        end
        R(a,k) = mean(tmpRefl); % average across band
    end
end
end

function X = maxminNDSI(V,S) %#ok<DEFNU>
% maximum (X(1)) and minimum (X(2)) values of NDSI for combinations of Visible
% and SWIR inputs

[VV,SS] = meshgrid(V,S);
N = (VV-SS)./(VV+SS);
X = [max(N(:)) min(N(:))];
end

function [NDSI,bands] = findBestNDSI(snowR,vnames,snowV,cloudV)

if ~isempty(snowR.cirrusBand) && any(strcmpi(vnames(cloudV),'cirrus'))
    sb = snowR.cirrusBand;
else
    sb = snowR.swirBands;
end
vb = snowR.visBands;
[vb,sb] = meshgrid(vb,sb);

[v,c] = meshgrid(snowV,cloudV);
vnameCombo = [v(:) c(:)];

% go through all cosines
WT = snowR.wholeTable;
cosIndex = WT.cosIndex;
vband = zeros(height(WT),size(vnameCombo,1));
sband = zeros(size(vband));
NDSI = table(cosIndex);
for h=1:height(WT)
    for combo=1:size(vnameCombo,1)
        sname = vnames{vnameCombo(combo,1)};
        cname = vnames{vnameCombo(combo,2)};
        tempT = WT(h,:);
        snowNDSI = (tempT.(sname)(vb)-tempT.(sname)(sb))./...
            (tempT.(sname)(vb)+tempT.(sname)(sb));
        cloudNDSI = (tempT.(cname)(vb)-tempT.(cname)(sb))./...
            (tempT.(cname)(vb)+tempT.(cname)(sb));
        Ndiff = snowNDSI-cloudNDSI;
        [row,col] = find(Ndiff==max(Ndiff(:)));
        vband(h,combo) = vb(row,col);
        sband(h,combo) = sb(row,col);
    end
end
vband = mode(vband(:));
sband = mode(sband(:));
for v=1:numel(snowV)
    vn = vnames{snowV(v)};
    xvb = WT{:,vn};
    NDSI.(vn) = (xvb(:,vband)-xvb(:,sband))./(xvb(:,vband)+xvb(:,sband));
end
for v=1:numel(cloudV)
    vn = vnames{cloudV(v)};
    xvb = WT{:,vn};
    NDSI.(vn) = (xvb(:,vband)-xvb(:,sband))./(xvb(:,vband)+xvb(:,sband));
end
bands = [vband sband];
end