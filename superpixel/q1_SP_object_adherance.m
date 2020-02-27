%% how well do superpixels adhear to the boundaires of known objects?
%Do some algorithms and paramaterizations outperform others?

%% LOAD DATA
%grab coincidence MODIS and Landsat aquisition
truthFilename = 'C:\raid\data\tstillinger\snowCloud\singleSceneTest\truthData\LC81480352013195LGN00_32_mask.png';
LS8C1_folder = 'C:\raid\data\tstillinger\snowCloud\singleSceneTest\LC08_L1TP_148035_20130714_20170503_01_T1';
MODISfilename = 'C:\raid\data\tstillinger\snowCloud\singleSceneTest\MOD09GA\MOD09GA.A2013195.h24v05.006.2015262122809.hdf';

tic %NEED To reproject modis data with same spatial footprint as originial data (500m foot print for cloud and bands)
[fullMODIS, subsetMODIS, fullOLI, subsetOLI, truth] = getTest_MODIS_and_OLI(truthFilename,MODISfilename,LS8C1_folder);
toc

%% SP scene reduction



%% BOUNDARY RECALL STATISTICS ON SUPERPIXEL CLOUD OBJECT ADHERENCE & HOW WELL DO SUPERPIXELS SEGMENT CLOUD FROM SNOW
% BASED ON: SP SIZE | BANDS USED TO SEGMENT (OR FEATURES) | ALGORITHM
%need to decide if you want to take boundaries of all superpixles and how
%well they follow cloud boundaries, or the singular boundary of all the
%superpixels that are majority cloud, as a boundary?

d =2;
classFlag=true;
X= subsetOLI.bands;
X(:,:,[8 9 10])=[]; %optical bands only, no cirrus
NDSIbands = [3 6];
[NDSI,pixNorm,distSAM,distE] = makeSPfeatures(X,NDSIbands);

%featureSets
set1=X(:,:,[6 5 3]);
set2=X(:,:,[6 5 4]);
set3=X(:,:,[4 3 2]);
set4= cat(3,NDSI,pixNorm,distSAM);
set5=cat(3,NDSI,distSAM,distE);

featureSets = cat(4,set1,set2,set3,set4,set5);
fSets = {'6 5 3','6 5 4','4 3 2','NDSI norm SAM','NDSI SAM Euc','SLICSAM'};
spSize = [10 100 500 1000 2500 5000 10000 100000]';
numRuns = size(featureSets,4)+1; %last run is SLICSAM

for k = 1:numRuns
    fset = repmat(fSets(k),size(spSize));
    boundaryRecallCloud=zeros(size(spSize));
    boundaryRecallSnow=zeros(size(spSize));
    UEcloud=zeros(size(spSize));
    UEsnow=zeros(size(spSize));
    
    if k<numRuns
        thisSet = featureSets(:,:,:,k);
    elseif k == numRuns %SLICSAM RUN - SLICSAM single set
        thisSet = X;
        m = .5; %m = k/100;%seems that m being 100x smaller thansize of super pixels might be a good size for it.
        SAMconst = 200; % this needs to be investigated (maybe) then pushed into function.
        nItr = 5;
    end
    
    
    for i = 1:length(spSize)
        [nR, nC,~] = size(thisSet);
        N = round((nR.*nC)./spSize(i));
        fillPxls = isnan(thisSet);
        validPixels = ~fillPxls;
        validPixels = all(validPixels,3);
        thisSet(fillPxls)=-9999;
        
        if k<numRuns
            [L,~] = superpixels(thisSet,N,'compactness',1,'method','slic');
            L(~validPixels)=0;
        elseif k == numRuns
            L = SAMSuperPixels(spSize(i),thisSet,m,SAMconst,nItr,validPixels);
        end
        
        [boundaryRecallCloud(i),classLcloud] = boundaryRecall(L,d,truth.cloud,classFlag);
        [boundaryRecallSnow(i),classLsnow] = boundaryRecall(L,d,truth.snow,classFlag);
            
        UEcloud(i) = underSegError(L,truth.cloud);
        UEsnow(i) = underSegError(L,truth.snow);
        
    end
    
    if k == 1
        fTbl = table(fset,spSize,boundaryRecallCloud,boundaryRecallSnow,UEcloud,UEsnow);
    else       
        tbl = table(fset,spSize,boundaryRecallCloud,boundaryRecallSnow,UEcloud,UEsnow);
        fTbl = [fTbl; tbl];
    end
end


%% plot results
spNames = unique(fTbl.fset);
figure;
hold on;
for k=1:length(spNames)
    rowIdcs = strcmpi(fTbl.fset,spNames{k});
    sizeSP = fTbl.spSize(rowIdcs);
    thisStat = fTbl.boundaryRecallCloud(rowIdcs);
    clr = rand(1,3);
    %plot((1000*1000)./sizeSP,BRclouds,'-','color',clr);%,sizeSP,BRsnow,'- -','color',clr);
    plot(sizeSP,thisStat,'-','color',clr);%,sizeSP,BRsnow,'- -','color',clr);
end
legend(spNames);
title('BR Cloud') 

figure;
hold on;
for k=1:length(spNames)
    rowIdcs = strcmpi(fTbl.fset,spNames{k});
    sizeSP = fTbl.spSize(rowIdcs);
    thisStat = fTbl.boundaryRecallSnow(rowIdcs);
    clr = rand(1,3);
    %plot((1000*1000)./sizeSP,BRclouds,'-','color',clr);%,sizeSP,BRsnow,'- -','color',clr);
    plot(sizeSP,thisStat,'-','color',clr);%,sizeSP,BRsnow,'- -','color',clr);
end
legend(spNames);
title('BR snow')

figure;
hold on;
for k=1:length(spNames)
    rowIdcs = strcmpi(fTbl.fset,spNames{k});
    sizeSP = fTbl.spSize(rowIdcs);
    thisStat = fTbl.UEcloud(rowIdcs);
    clr = rand(1,3);
    %plot((1000*1000)./sizeSP,BRclouds,'-','color',clr);%,sizeSP,BRsnow,'- -','color',clr);
    plot(sizeSP,thisStat,'-','color',clr);%,sizeSP,BRsnow,'- -','color',clr);
end
legend(spNames);
title('UE Cloud')

figure;
hold on;
for k=1:length(spNames)
    rowIdcs = strcmpi(fTbl.fset,spNames{k});
    sizeSP = fTbl.spSize(rowIdcs);
    thisStat = fTbl.UEsnow(rowIdcs);
    clr = rand(1,3);
    %plot((1000*1000)./sizeSP,BRclouds,'-','color',clr);%,sizeSP,BRsnow,'- -','color',clr);
    plot(sizeSP,thisStat,'-','color',clr);%,sizeSP,BRsnow,'- -','color',clr);
end
legend(spNames);
title('UE Snow')


%% STATISTICS ON SUPERPIXEL CLOUD OBJECT ADHERENCE & HOW WELL DO SUPERPIXELS SEGMENT CLOUD FROM SNOW
% BASED ON: SP SIZE | BANDS USED TO SEGMENT (OR FEATURES) | ALGORITHM

X= subsetOLI.bands;
X(:,:,[8 9 10])=[]; %optical bands only, no cirrus
NDSIbands = [3 6];
[NDSI,pixNorm,distSAM,distE] = makeSPfeatures(X,NDSIbands);

%featureSets
set1=X(:,:,[6 5 3]);
set2=X(:,:,[6 5 4]);
set3=X(:,:,[4 3 2]);
set4= cat(3,NDSI,pixNorm,distSAM);
set5=cat(3,NDSI,distSAM,distE);

featureSets = cat(4,set1,set2,set3,set4,set5);
fSets = {'6 5 3','6 5 4','4 3 2','NDSI norm SAM','NDSI SAM Euc'};
spSize = [100 500 1000 2500 5000 10000 100000]';

for k = 1:size(featureSets,4)
    fset = repmat(fSets(k),size(spSize));
    thisSet = featureSets(:,:,:,k);
    iReduceClouds=zeros(size(spSize));
    RMSE=zeros(size(spSize));
    MAE=zeros(size(spSize));
    bias=zeros(size(spSize));
    fMixedClouds=zeros(size(spSize));
    medianMixClouds=zeros(size(spSize));
    medianaAreaClouds=zeros(size(spSize));
    iReduceSnow=zeros(size(spSize));
    fMixedSnow=zeros(size(spSize));
    medianMixSnow=zeros(size(spSize));
    medianAreaSnow=zeros(size(spSize));
    
    for i = 1:length(spSize)
        [nR, nC,~] = size(thisSet);
        N = round((nR.*nC)./spSize(i));
        fillPxls = isnan(thisSet);
        validPixels = ~fillPxls;
        validPixels = all(validPixels,3);
        thisSet(fillPxls)=-9999;
        
        %try with SLIC and SLICSAM
        tic
        [L,~] = superpixels(thisSet,N,'compactness',1,'method','slic');
        toc
        L(~validPixels)=0;
        
        %fractional cover of truth classes within each superpixel
        [featureSet] = calcSuperPixelFeatures(L,truth.snow,truth.cloud,truth.neither);
        cloudFrac = featureSet(:,:,2);
        snowFrac = featureSet(:,:,1);
        F1=cloudFrac;
        T1=truth.cloud;
        F2=snowFrac;
        T2=truth.snow;
        [iReduceClouds(i),RMSE(i),MAE(i),bias(i),fMixedClouds(i),medianMixClouds(i),medianaAreaClouds(i),iReduceSnow(i),fMixedSnow(i),medianMixSnow(i),medianAreaSnow(i)] = superPixelShapeStats(L,F1,T1,F2,T2);
        
    end
    
    if k == 1
        fTbl = table(fset,spSize,iReduceClouds,RMSE,MAE,bias,fMixedClouds,medianMixClouds,iReduceSnow,fMixedSnow,medianMixSnow);
        figure
        plot(tbl.spSize,tbl.MAE)
        hold on
    else
        tbl = table(fset,spSize,iReduceClouds,RMSE,MAE,bias,fMixedClouds,medianMixClouds,iReduceSnow,fMixedSnow,medianMixSnow);
        plot(tbl.spSize,tbl.MAE)
        fTbl = [fTbl; tbl];
    end
end
hold off
legend('6 5 3','6 5 4','4 3 2','NDSI norm SAM','NDSI SAM Euc')
xlim([0 10000])




%% boundary recall & undersegmentation error - counting cloud shadows as clearsky for now
%BOUNDARY RECALL _ NEED A WAY TO NOT PENALIZE FOR BOUNDARIES WITHIN THE CLOUD
%For these, focusing on accurate segmentation of clouds in the image.
%find a way to limit boundary recall to just the

cloudTruth = truth.cloud;
clearSkyTruth = ~truth.cloud;

%only calculate statistics on superpixels that touch clouds?
d =2;
classFlag = true;
[recall,classL] = boundaryRecall(L,d,truth.cloud,classFlag);

%% BOUNDARY RECALL FIGURE FOR PRESENTATION
fRGB = landsat8fRGB(X, [6,5,3]);


%truth boundary and superpixel boundary
truthBM = boundarymask(truth.cloud,8);
SE = strel('square',d);
truthBM = imdilate(truthBM,SE);
truthBoundary = cat(3, zeros(size(L)), zeros(size(L)), zeros(size(L)));

SPB = boundarymask(classLcloud,4);
SPBoundary = cat(3, ones(size(L)), zeros(size(L)), zeros(size(L)));

figure
imshow(fRGB)
hold on
h = imshow(truthBoundary);
H = imshow(SPBoundary);
hold off
set(h, 'AlphaData', truthBM)
set(H, 'AlphaData', SPB)

figure
%%




figure
imshow(fRGB)
hold on
h = imshow(blkSPbounary);
hold off
set(h, 'AlphaData', BW)

%over truth mask - good for talk
blkSPbounary = cat(3, ones(size(L)), zeros(size(L)), zeros(size(L)));
figure
imshow(truth.cloud)
hold on
h = imshow(blkSPbounary);
hold off
set(h, 'AlphaData', BW)


[UE] = underSegError(L,cloudTruth,clearSkyTruth);

%%
% %
% %
% % figure
% % histogram(snowFrac(cloudAndSnow),'Normalization','probability',5)
% % hold on
% % histogram(cloudFrac(cloudAndSnow),'Normalization','probability',5)
% % hold off
% % legend('snow','cloud')
% %
% %
% % snowFrac = featureSet(:,:,1);
% % snow = single(snowFrac(idxS(iS));
% %
% %
% %
% % snow = single(snowFrac(idxS));
% %
% % neitherFrac = featureSet(:,:,3);
% % figure; histogram(cloudFrac(cloudSPs),'Normalization','probability')
% % mixups = cloudFrac>40&snowFrac>40;
% % figure; histogram(cloudFrac(cloudFrac>0))
% % figure; histogram(snowFrac(snowFrac>0))
% %
% % figure; histogram(cloudFrac(mixups))
% % figure; histogram(snowFrac(mixups))
% %
% %
% %
% % mixupsMask = labeloverlay(fRGB,mixups,'IncludedLabels',1,'colormap',[0.15,0.15,0.15],'Transparency',0.25);
% % figure;imshow(mixupsMask)

%% SPECTRAL VARIABILITY OF SUPERPIXELS
[C] = calcSuperPixelMeans(L,X);



%% VISUALIZATIONS
fRGB = landsat8fRGB(X, [6,5,3]);

BW = boundarymask(L,4);
blkSPbounary = cat(3, zeros(size(L)), zeros(size(L)), zeros(size(L)));
figure
imshow(fRGB_M)
hold on
h = imshow(blkSPbounary);
hold off
set(h, 'AlphaData', BW)
