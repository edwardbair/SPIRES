% How do various methods of generating superpixels impact there ability to
% segment snow and clouds?
% 
% picking one exmample to test on - 1000x1000 usgs subscene with truth classes.
%looking at MODIS and Landsat - landsat data is TOA, modis data is Rs for
%now.

%grab coincidence MODIS and Landsat aquisition
truthFilename = 'C:\raid\data\tstillinger\snowCloud\singleSceneTest\truthData\LC81480352013195LGN00_32_mask.png';
LS8C1_folder = 'C:\raid\data\tstillinger\snowCloud\singleSceneTest\LC08_L1TP_148035_20130714_20170503_01_T1';
MODISfilename = 'C:\raid\data\tstillinger\snowCloud\singleSceneTest\MOD09GA\MOD09GA.A2013195.h24v05.006.2015262122809.hdf';

tic %NEED To reproject modis data with same spatial footprint as originial data (500m foot print for cloud and bands)
[fullMODIS, subsetMODIS, fullOLI, subsetOLI, truth] = getTest_MODIS_and_OLI(truthFilename,MODISfilename,LS8C1_folder);
toc

%% Visualize operational cloud masks

%false color images
fRGB_500m = rgbMODIS( fullMODIS.bands,[6 2 1] );
fRGB_30m = landsat8fRGB(fullOLI.bands, [6,5,3]);
fRGB_500m_sub = rgbMODIS( subsetMODIS.bands,[6 2 1] );
fRGB_30m_sub = landsat8fRGB(subsetOLI.bands, [6,5,3]);


figure
subplot(2,2,1)
imshow(fRGB_30m)
title('full scene OLI fRGB')
subplot(2,2,2)
imshow(fRGB_500m)
title('full OLI scene MODIS')
subplot(2,2,3)
imshow(fRGB_30m_sub)
title('1000x1000 oli subscene')
subplot(2,2,4)
imshow(fRGB_500m_sub)
title('MODIS subset to 1000x1000 oli subscene')


%cloud masks over the false  color images
modisCloud = fullMODIS.cloudMask == 1 ;%| fullMODIS.cloudMask == 2 | fullMODIS.cirrusMask == 2 |  fullMODIS.cirrusMask == 3;
modisCloudMask = labeloverlay(fRGB_500m,modisCloud,'IncludedLabels',1,'colormap',[0.15,0.15,0.15],'Transparency',0.25);

modisCloud_sub = subsetMODIS.cloudMask == 1 ;%| subsetMODIS.cloudMask == 2| subsetMODIS.cirrusMask == 2 |  subsetMODIS.cirrusMask == 3;
modisCloudMask_sub = labeloverlay(fRGB_500m_sub,modisCloud_sub,'IncludedLabels',1,'colormap',[0.15,0.15,0.15],'Transparency',0.25);

oliCloud = fullOLI.cloudMask; %== 1 | fullMODIS.cloudMask == 2;
oliCloudMask = labeloverlay(fRGB_30m,oliCloud,'IncludedLabels',1,'colormap',[0.15,0.15,0.15],'Transparency',0.25);

oliCloud_sub = subsetOLI.cloudMask;% == 1 | subsetMODIS.cloudMask == 2;
oliCloudMask_sub = labeloverlay(fRGB_30m_sub,oliCloud_sub,'IncludedLabels',1,'colormap',[0.15,0.15,0.15],'Transparency',0.25);


figure
subplot(2,2,1)
imshow(oliCloudMask)
title('full scene OLI fRGB')
subplot(2,2,2)
imshow(modisCloudMask)
title('full OLI scene MODIS')
subplot(2,2,3)
imshow(oliCloudMask_sub)
title('1000x1000 oli subscene')
subplot(2,2,4)
imshow(modisCloudMask_sub)
title('MODIS subset to 1000x1000 oli subscene')

%% visualize ned's current cloud mask
%cloud QAbits == 1 & band6>0.2
spireCloudV1 = modisCloud & fullMODIS.bands(:,:,6)>0.2;
spireCloudV1_sub = modisCloud_sub & subsetMODIS.bands(:,:,6)>0.2;
spireCloudMaskv1 = labeloverlay(fRGB_500m,spireCloudV1,'IncludedLabels',1,'colormap',[0.15,0.15,0.15],'Transparency',0.25);
spireCloudMaskv1_sub = labeloverlay(fRGB_500m_sub,spireCloudV1_sub,'IncludedLabels',1,'colormap',[0.15,0.15,0.15],'Transparency',0.25);

figure
subplot(2,2,1)
imshow(oliCloudMask)
title('full scene OLI fRGB')
subplot(2,2,2)
imshow(spireCloudMaskv1)
title('full OLI scene MODIS w/ spirev1 cloudMask')
subplot(2,2,3)
imshow(oliCloudMask_sub)
title('1000x1000 oli subscene')
subplot(2,2,4)
imshow(spireCloudMaskv1_sub)
title('MODIS subset to 1000x1000 oli subscene w/ spirev1 cloudMask')

%% run various superpixel segmentation algorihtms.


% NDSI | norm | spectral angle

%       vs

% SWIR band | NIR band | Vis Band (frgb)


%        with

%    SLIC  or   SLICSAM


%see which one does the best job of putting nsow and lcoud pixels in
%different superpixels. 

%% features other than band reflectance for use in superpixel segmentation
xL= subsetOLI.bands;
xL(:,:,8:10)=[];%not using thermal data or cirrus band for now (same bands as modis)
NDSIbands = [3 6];

xM= fullMODIS.bands;
NDSIbands = [4 6];


[NDSI,pixNorm,distSAM,distE] = makeSPfeatures(xL,NDSIbands);


%visualize
XL = cat(3,NDSI,distSAM,distE);
fRGB = landsat8fRGB(XL, [1,2,3]);
%% superpixel size
spSize = [100 500 1000 2500 5000 10000 100000];
 %500m modis pixel = ~16.6x16.6 landsat pixels = 275.56 pixels (~280)
 %250m modis pixel = ~8.3x8.3 landsat pixels = 68.89 pixels    (~70)

%% built in MATLAB superpixels (Normal SLIC alg)

superPixelSize = 250; %MAKE SURE THIS IS CORRECT
[nR, nC,~] = size(XL);
N = round((nR.*nC)./superPixelSize);

fillPxls = isnan(XL);
validPixels = ~fillPxls;
% % % %erode by 1 to avoid texture flags along boarder. - alterante fix could be to transpose image around its boarder.
% % % SE=strel('square',8); 
% % % validPixels = imerode(validPixels,SE);
validPixels = all(validPixels,3);
% % % 
% % % figure; imshow(validPixels)
XL(fillPxls)=-9999;
% % % A = X(:,:,[6 5 3]);
% % % %convert colorspace to lab
%lab = rgb2lab(A);


%try with SLIC and SLICSAM
A =XL;
tic
[L,~] = superpixels(A,N,'compactness',1,'method','slic');
toc
L(~validPixels)=0;
 
BW = boundarymask(L,4);
blkSPbounary = cat(3, zeros(size(L)), zeros(size(L)), zeros(size(L)));
figure
imshow(fRGB_30m_sub)
hold on
h = imshow(blkSPbounary);
hold off
set(h, 'AlphaData', BW)


%% FEATURES OF THE SUPERPIXELS
%NDSI,pixNorm,distSAM,distE;
 [featureSet] = calcSuperPixelFeatures(L,truth.snow,truth.cloud,truth.neither);
 %[C] = calcSuperPixelMeans(L,X);

 cloudFrac = featureSet(:,:,2);
 snowFrac = featureSet(:,:,1);
 neitherFrac = featureSet(:,:,3);
 figure; histogram(cloudFrac)
 
 mixups = cloudFrac>0.4&snowFrac>0.4;
 figure; histogram(cloudFrac(cloudFrac>0))
figure; histogram(snowFrac(snowFrac>0))

figure; histogram(cloudFrac(mixups))
figure; histogram(snowFrac(mixups))



mixupsMask = labeloverlay(fRGB_30m_sub,mixups,'IncludedLabels',1,'colormap',[0.15,0.15,0.15],'Transparency',0.25);
figure;imshow(mixupsMask)

%% 3 bands - RGB (2,3,4)

%3 bands fRGB(6, 5, 3)

%NDSI | Norm | specFlatness



%my multiband SLICSAM

%all bands

%all bands + norm

%allbands + NDSI | Norm | specFlatness




%figure showing all




%compare output to truth mask

% average superpixel % cover of cover class of interest (cloud or snow)
%(superpixels that include at
% least 1 cloud mask pixel)

