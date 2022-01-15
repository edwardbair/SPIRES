function [likelySnow, maybeSnow, likelyCloud, maybeCloud, varargout ] =...
    filterCloudSnow(snowR, surfaceRefl, solarZ)
% [likelySnow, maybeSnow, likelyCloud, maybeCloud, [maybeCirrus] ] = filterCloudSnow(snowR, surfaceRefl, solarZ)
%
% General function for snow/cloud filtering in multispectral data, and used
% by the MODIS and Landsat cloud filters
%
% Input
%   snowR - snow reflectivity structure from prepCloudSnowFilter, which
%       contains values for coarse snow, fine snow, dirty snow, thick and
%       thin clouds, and cirrus
%   surfaceRefl - reflectance product (for example MOD09) in the same
%       bands over which snowR has been calculated, e.g. MxNx7 for MODIS
%   solarZ - solar zenith angles covering the same image, 2D (could be
%       scalar or different 2D grid size than surfaceRefl)
%
% Output
%   likelySnow - bright enough to be snow, but not bright enough in SWIR
%       for cloud
%   maybeSnow - bright enough to be dirty snow, but not a cloud
%   likelyCloud - bright enough to be snow or thick cloud, but too bright
%       in SWIR to be snow
%   maybeCloud - probably a thin cloud
% Optional output
%   maybeCirrus - but only if the snowR structure includes a cirrus band

numArg = 4;
if nargout>numArg && ~isempty(snowR.cirrusBand)
    doCirrus = true;
else
    doCirrus = false;
end

cosZ = cosd(solarZ);
cosLimits = [nanmin(cosZ(:)) nanmax(cosZ(:))];
% fill NaNs if any
if any(isnan(cosZ(:)))
    cosZ = inpaint_nans(double(cosZ),4);
end
% resize if needed
imageSize = [size(surfaceRefl,1) size(surfaceRefl,2)];
if isscalar(cosZ)
    cosZ = repmat(cosZ,imageSize);
elseif ~isequal(size(cosZ),imageSize)
    cosZ = imresize(cosZ,imageSize);
end
cosZ = truncateLimits(cosZ,cosLimits);

% multiply cosines by multiplier, round to integers
cosZ = round(snowR.cosineMultiplier*cosZ);
cosZ(cosZ<=0) = 1; % eliminate points at or below horizon
cosZ(cosZ>snowR.cosineMultiplier) = snowR.cosineMultiplier; % shouldn't happen, but just in case

% propagate model values in the snowR structure to the image, based on
% actual cosZ values
if doCirrus
    tableVariable = {'fineSnow','coarseSnow','dirtySnow','thickCloud','thinCloud','cirrus'};
else
    tableVariable = {'fineSnow','coarseSnow','dirtySnow','thickCloud','thinCloud'};
end
modelAllBands = propagateModel(tableVariable,cosZ,snowR.wholeTable,class(surfaceRefl));
modelNDSI = propagateModel(tableVariable,cosZ,snowR.NDSI,class(surfaceRefl));

vb = snowR.ndsiBands(1);
sb = snowR.ndsiBands(2);
actualNDSI = (surfaceRefl(:,:,vb)-surfaceRefl(:,:,sb))./...
    (surfaceRefl(:,:,vb)+surfaceRefl(:,:,sb));
if doCirrus
    cb = snowR.ndsiBands(3);
    cirrusNDSI = (surfaceRefl(:,:,vb)-surfaceRefl(:,:,cb))./...
        (surfaceRefl(:,:,vb)+surfaceRefl(:,:,cb));
end

% some measures for snow, cloud, cirrus need to be pre-allocated
maskVariable = {'snowVis','snowNIR','dirtySnowVis','dirtySnowNIR',...
    'cloudVis','cloudNIR','cloudSWIR','tcVis','tcNIR','tcSWIR'};
for v=1:length(maskVariable)
    mask.(maskVariable{v}) = false(imageSize);
end

% any visible band bright enough
for b=snowR.visBands
    thisRefl = surfaceRefl(:,:,b);
    mask.snowVis = mask.snowVis |...
        (thisRefl>=modelAllBands.coarseSnow(:,:,b));
    mask.dirtySnowVis = mask.dirtySnowVis |...
        (thisRefl>=modelAllBands.dirtySnow(:,:,b));
    mask.cloudVis = mask.cloudVis |...
        (thisRefl>=modelAllBands.thickCloud(:,:,b));
    mask.tcVis = mask.tcVis |...
        (thisRefl>=modelAllBands.thinCloud(:,:,b));
end
% any NIR band
for b=snowR.nirBands
    thisRefl = surfaceRefl(:,:,b);
    mask.snowNIR = mask.snowNIR |...
        (thisRefl>=modelAllBands.coarseSnow(:,:,b));
    mask.dirtySnowNIR = mask.dirtySnowNIR |...
        (thisRefl>=modelAllBands.dirtySnow(:,:,b));
    mask.cloudNIR = mask.cloudNIR |...
        (thisRefl>=modelAllBands.thickCloud(:,:,b));
    mask.tcNIR = mask.tcNIR |...
        (thisRefl>=modelAllBands.thinCloud(:,:,b));
end
% any SWIR band
for b=snowR.swirBands
    thisRefl = surfaceRefl(:,:,b);
    mask.cloudSWIR = mask.cloudSWIR |...
        (thisRefl>=modelAllBands.thickCloud(:,:,b));
    mask.tcSWIR = mask.tcSWIR |...
        (thisRefl>=modelAllBands.thinCloud(:,:,b));
end
% NDSI
mask.snowNDSI = actualNDSI>=modelNDSI.fineSnow;
mask.cloudNDSI = actualNDSI<=modelNDSI.thickCloud &...
    actualNDSI>=modelNDSI.thinCloud;
mask.tcNDSI = actualNDSI<modelNDSI.thinCloud & actualNDSI>0;
if doCirrus
    mask.cirrusNDSI = cirrusNDSI<=modelNDSI.cirrus;
end

% assign likelihoods
cloudScore = mask.cloudVis+mask.cloudNIR+2*mask.cloudSWIR+mask.cloudNDSI+...
    mask.tcVis+mask.tcNIR+mask.tcSWIR+mask.tcNDSI-mask.snowNDSI;
likelyCloud = cloudScore>2;
maybeCloud = cloudScore>0 & ~likelyCloud;
if doCirrus
    cScore = mask.cirrusNDSI +...
        (surfaceRefl(:,:,cb)>=modelAllBands.cirrus(:,:,cb));
    maybeCirrus = cScore>0;
    varargout{1} = maybeCirrus;
end
snowScore = mask.snowVis+mask.snowNIR+mask.dirtySnowVis+mask.dirtySnowNIR+mask.snowNDSI;
likelySnow = snowScore>3;
maybeSnow = snowScore>2 & ~likelySnow;
% fix pixels identified as both likely snow and cloud
t = likelySnow & likelyCloud;
gs = snowScore>cloudScore;
likelyCloud(t&gs) = false;
maybeCloud(t&gs) = true;
likelySnow(t & ~gs) = false;
maybeSnow(t & ~gs) = true;

end

function S=propagateModel(tableVariable,cosZ,whichTable,classname)

% size of variables in output structure
vSiz = [size(cosZ) length(whichTable{1,tableVariable{1}})];
U = unique(cosZ(:));

% reshape for faster indexing
cosZ = reshape(cosZ,numel(cosZ),1);
for k=1:length(tableVariable)
    % reshape cube to support faster indexing
    numPix = vSiz(1)*vSiz(2);
    X = zeros(numPix,vSiz(3),classname);
    % values in array for indexing in parfor loop
    values = zeros(length(U),length(whichTable{1,tableVariable{k}}));
    for c=1:length(U)
        v = whichTable{:,'cosIndex'}==U(c);
        values(c,:) = whichTable{v,tableVariable{k}};
    end
    for band=1:size(X,2)
        Y = zeros(numPix,1,classname);
        for c=1:length(U)
            t = cosZ==U(c);
            Y(t) = values(c,band);
        end
        X(:,band) = Y;
    end
    S.(tableVariable{k}) = reshape(X,vSiz(1),vSiz(2),vSiz(3));
end
end