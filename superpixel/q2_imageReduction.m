truthFilename = 'C:\raid\data\tstillinger\snowCloud\singleSceneTest\truthData\LC81480352013195LGN00_32_mask.png';
LS8C1_folder = 'C:\raid\data\tstillinger\snowCloud\singleSceneTest\LC08_L1TP_148035_20130714_20170503_01_T1';
MODISfilename = 'C:\raid\data\tstillinger\snowCloud\singleSceneTest\MOD09GA\MOD09GA.A2013195.h24v05.006.2015262122809.hdf';

%run SLICSAM
%uniquetol 5% - sompare image size before/after
%compare image size snow/cloud images
%
runMODISFlag=false;
modelMODISFlag=false;
[fullMODIS, subsetMODIS, fullOLI, subsetOLI, truth] = getTest_MODIS_and_OLI(truthFilename,MODISfilename,LS8C1_folder,runMODISFlag,modelMODISFlag);


%%

X= subsetOLI.bands;
cloudTruth = truth.cloud;
snowTruth = truth.snow;
spSize = [10 100 500 1000 2500 5000 10000 100000]';

X(:,:,[8 9 10])=[]; %optical bands only, no cirrus

thisSet = X;
m = .5; %m = k/100;%seems that m being 100x smaller thansize of super pixels might be a good size for it.
SAMconst = 200; % this needs to be investigated (maybe) then pushed into function.
nItr = 5;

for aa = 1:length(spSize)
    [nR, nC,~] = size(thisSet);
    N = round((nR.*nC)./spSize(aa));
    fillPxls = isnan(thisSet);
    validPixels = ~fillPxls;
    validPixels = all(validPixels,3);
    thisSet(fillPxls)=-9999;
    L = SAMSuperPixels(spSize(aa),thisSet,m,SAMconst,nItr,validPixels);
    L(~validPixels)=0;
    fillIdx = L ==-1;
    L(fillIdx)=0;%HACK
end


%figure of the example
BW = boundarymask(L,4);
fRGB = landsat8fRGB(X, [6,5,3]);
blkSPbounary = cat(3, zeros(size(L)), zeros(size(L)), zeros(size(L)));
figure
imshow(fRGB)
hold on
h = imshow(blkSPbounary);
hold off
set(h, 'AlphaData', BW)

%rearagne superpixel mean values to 2D
[spMeans] = calcSuperPixelMeans(L,X);
spMeans = round(spMeans,2);
%reshape all pixles to 2D values
[r,c,d]=size(X);
pxlVals = reshape(X,[r*c,d]);
pxlVals = round(pxlVals,2);

%compare histograms of pixel values

tol = 0.05;
utSP = uniquetol(spMeans,tol,'ByRows',true);
 

utPxl = uniquetol(pxlVals,tol,'ByRows',true);
whos spMeans C utPxl pxlVals utSP

%full scene redution
whos pxlVals spMeans utSP

%SNOW Reduction
snowSPMeans=spMeans;

deleteSPs=unique(L(~truth.snow));
snowSPMeans(deleteSPs,:)=[];
utSnowSP = uniquetol(snowSPMeans,tol,'ByRows',true);
%reshape all pixles to 2D values
snowpxlVals = pxlVals;
snowpxlVals(~truth.snow(:),:)=[];
whos snowpxlVals snowSPMeans utSnowSP 



%FIG 1 - normalized histogram
figure
xText = {'band 1','band 2','band 3','band 4','band 5','band 6','band 7'};
for i = 1:7
subplot(2,4,i)
histogram(snowpxlVals(:,i),10,'Normalization','probability','FaceAlpha',.5)
hold on
histogram(snowSPMeans(:,i),10,'Normalization','probability','FaceAlpha',.5)
histogram(utSnowSP(:,i),10,'Normalization','probability','FaceAlpha',.5)
xlabel(xText{i})
set(gca,'fontsize',15)
xlim([0 1.2])
end
legend('all snow pixels','all snow superpixels','unique snow superpixels')

figure
xText = {'band 1','band 2','band 3','band 4','band 5','band 6','band 7'};
for i = 1:7
subplot(2,4,i)
histogram(utSnowSP(:,i),'Normalization','probability','FaceAlpha',.5)
hold on
histogram(snowSPMeans(:,i),'Normalization','probability','FaceAlpha',.5)
histogram(snowpxlVals(:,i),'Normalization','probability','FaceAlpha',.5)
xlabel(xText{i})
set(gca,'fontsize',15)
xlim([0 1.25])
end
legend('unique snow superpixels','all snow superpixels','all snow pixels')


% % 
% % if cntr == 1
% %     fTbl = table();
% % else
% %     tbl = table();
% %     fTbl = [fTbl; tbl];
% % end
% % cntr = cntr+1;


