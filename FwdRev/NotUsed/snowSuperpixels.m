function S = snowSuperpixels(wavelength,reflectance,elevation,cosZ,muS,Ftrans,varargin)
% S = snowSuperpixels(wavelength,reflectance,elevation,cosZ,muS,Ftrans [,viewF])
%create superpixels for hyperspectral snow scene
%
%Input
% wavelength - wavelength, assume nm
% reflectance - image reflectance that has been through correctAVNG in BIP
%   format (use bsq2bip to create if necessary)
% elevation - elevation grid same size as 2nd and 3rd reflectance dimensions
% cosZ - scalar or grid, cosine of solar zenith
% muS - grid of cosines of illumination on slopes (can be corrected for
%   horizons by setting to zero where H>asind(cosZ))
% Ftrans - atmosphere transmittance function, squared
%Optional input
% viewF - view factor (not used to find superpixels, but keeps track for
%   subsequent calculations
%
%Output, structure S with the following fields
% spectralSP - mean spectrum of each superpixel
% spectralWeightedSP - mean weighted (squared atmospheric transmission)
%   spectrum of each pixel
% pixIdx - mapping of superpixel IDs to image
% L - superpixel labels
% nS - number of superpixels
% rgbImg - superpixel 3 band image, histogram equalized, with R = norm of
%   spectrum for each pixel, G = spectral angle between spectrum and clean
%   snow spectrum, B = cosine of illumination angle on slope
% superpixelImg - superpixels corresponding to rgbImg

p = inputParser;
addRequired(p,'wavelength',@(x) isvector(x) && isnumeric(x) && all(x>0,'all'))
addRequired(p,'reflectance',@(x) ndims(x)==3 && isnumeric(x) &&...
    all(x(~isnan(x))>=0,'all'))
addRequired(p,'elevation',@(x) isnumeric(x) && (isscalar(x) || ismatrix(x)))
addRequired(p,'cosZ',@(x) isnumeric(x) && all(x>0,'all') &&...
    all(x<=1,'all') && (isscalar(x) || ismatrix(x)))
addRequired(p,'muS',@(x) isnumeric(x) && all(x>=0,'all') &&...
    all(x<=1,'all') && ismatrix(x))
addRequired(p,'Ftran',@(x) strcmpi(class(x),'cfit') ||...
    contains(class(x),'interpolant','IgnoreCase',true))
addOptional(p,'viewF',1,@(x) isnumeric(x) && all(x>0,'all') &&...
    all(x<=1,'all') && (isscalar(x) || ismatrix(x)))
parse(p,wavelength,reflectance,elevation,cosZ,muS,Ftrans,varargin{:});

% atmospheric transmittance (squared)
wt = getAtmosProp(Ftrans,wavelength,mean(elevation,'all','omitnan'),mean(cosZ,'all','omitnan'));
assert(isequal(length(wavelength),size(reflectance,1)),...
    'number of wavelengths (%d) must be same as 1st dimension of reflectance (%d)',...
    length(wavelength),size(reflectance,1))
cubeWt = repmat(wt,1,size(reflectance,2),size(reflectance,3));

% reflectance of clean, fine snow
P = fwdPrescription('snow','cosZ',mean(cosZ,'all','omitnan'),'radius',100,'wavelength',...
    wavelength,'waveu','nm','elevation',mean(elevation,'all','omitnan'));
Rsnow =  SPIReS_fwd(P);

% scaled norm of reflectance
trefNaN = isnan(reflectance);
reflectance(trefNaN) = 0;
Nraw = squeeze(vecnorm(cubeWt.*reflectance));

% identify superpixels
snowCube = repmat(Rsnow,1,size(reflectance,2),size(reflectance,3));
cosA = squeeze(dot(cubeWt.*reflectance,cubeWt.*snowCube))./...
    (Nraw.*squeeze(vecnorm(cubeWt.*snowCube)));
t = isnan(cosA);
muS(t) = NaN;
rgbImg = uint16(round(double(intmax('uint16'))*cat(3,rescale(Nraw),cosA,muS)));
for b=1:3
    X = rgbImg(:,:,b);
    rgbImg(:,:,b) = histeq(X);
end
[L,nS] = superpixels(rgbImg,round(numel(rgbImg)/100));

% rbg bands in the superpixels
outputImg = zeros(size(rgbImg),'like',rgbImg);
idx = label2idx(L);
numRows = size(rgbImg,1);
numCols = size(rgbImg,2);
for labelVal = 1:nS
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    outputImg(redIdx) = round(mean(rgbImg(redIdx)));
    outputImg(greenIdx) = round(mean(rgbImg(greenIdx)));
    outputImg(blueIdx) = round(mean(rgbImg(blueIdx)));
end

% spectra and topography in the superpixels
spectralSP = zeros(length(wavelength),nS);
spectralWeightedSP = zeros(size(spectralSP));
viewF = p.Results.viewF;
doTopography = ~(isscalar(elevation) && isscalar(muS) && isscalar(cosZ) &&...
    isscalar(viewF));
if doTopography
    topography = zeros(5,nS); % 1 elevation, 2 cosZ, 3 muS, 4 viewF, 5 nPixel
end
reflectance(trefNaN) = NaN;
for labelVal = 1:nS
    pixIdx = idx{labelVal};
    spectralSP(:,labelVal) = mean(reflectance(:,pixIdx),2,'omitnan');
    spectralWeightedSP(:,labelVal) = mean(cubeWt(:,pixIdx).*reflectance(:,pixIdx),2,'omitnan');
    nPixel = numel(pixIdx);
    % topography
    if doTopography
        if isscalar(elevation)
            meanElev = elevation;
        else
            meanElev = mean(elevation(pixIdx),'omitnan');
        end
        if isscalar(cosZ)
            meanCosine = cosZ;
        else
            meanCosine = mean(cosZ(pixIdx),'omitnan');
        end
        if isscalar(muS)
            mean_muS = muS;
        else
            mean_muS = mean(muS(pixIdx),'omitnan');
        end
        if isscalar(viewF)
            meanView = viewF;
        else
            meanView = mean(viewF(pixIdx),'omitnan');
        end
        topography(:,labelVal) = [meanElev; meanCosine; mean_muS; meanView; nPixel];
    else
        topography(:,labelVal) = [elevation; cosZ; muS; viewF; nPixel];
    end
end

tnan = isnan(sum(spectralWeightedSP,1));
spectralSP(:,tnan) = [];
spectralWeightedSP(:,tnan) = [];
topography(:,tnan) = [];

% output into structure
S.wavelength = wavelength(:).';
S.spectralSP = spectralSP.';
S.spectralWeightedSP = spectralWeightedSP.';
S.idx = idx;
S.L = L;
S.rgbImg = rgbImg;
S.superpixelImg = outputImg;
S.topography = topography.';
S.topoID = {'elevation','cosSolarZ','cosIllumSlope','viewFactor','numPixels'};
end