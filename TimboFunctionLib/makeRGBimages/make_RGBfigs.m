%print two images in each plot
%Image 1
%        false color rgb of the scene the other
%Image 2
%       anything not snow - black
%       snow - turqoise
%       snow confused with clouds - dark blue
%       clouds - white
%       clouds confused with snow - grey
% subplot 3?
%       stats for the results -
%       unsupervised cluster results
%       fraction of pixelsin each class missclassifed
%       current cloud maks error statistics

%Location of all truth masks
'C:\raid\data\tstillinger\snowCloud\truthMasks_matlab';

%% test scene   LC80470252014226LGN00

%truthMask
load('C:\raid\data\tstillinger\snowCloud\truthMasks_matlab\LC80470252014226LGN00_USGS_truthMasks.mat');
%similarity data
load('C:\raid\scratch\timbo\snowCloud2017\similarityCalcs\LC80470252014226LGN00_5percentSimilar.mat');

%load('C:\raid\scratch\timbo\snowCloud2017\similarityCalcs\onePercent\USGS\LC80470252014226LGN00_1percentSimilar.mat');

%TOAdata
tiffFolder = 'C:\raid\data\tstillinger\snowCloud\snowCloudTruthMasks_Raw\ManualMasks_USGSanalyist\l8cloudmasks\clouds_snow_ice_masks';
geotiffFilename = 'LC80470252014226LGN00_11_data.tif';

[ Rtoa ] = USGSmaskTOA( fullfile(tiffFolder,geotiffFilename) );



% SNOWIDX and CLOUDIDX is the idex location of the PIXEL in the
% CloudPixSpec and SnowPixSepc vectors! NOT the locaiton of the pixel in
% the image....

[similarSnow, similarCloud] = similarMasks_USGS(snowIdx,cloudIdx,cloudPxlSpec,snowPxlSpec);

figure
imshow(similarSnow)
title('similar snow')
figure
imshow(similarCloud)
title('similar cloud')



% print image of scene



% color the truth masks

% color the 5% similar spectra in each image.


[ maskRGB, cmap, labels ] = maskRGBs( truthMask );
numLabels = length(labels);

figure
imshow(maskRGB,cmap)
h = colorbar;
h.Ticks = 0.5:1:numLabels-0.5;
h.TickLabels = labels';

%add catagories
% 1% AND 5% similar cloud and snow (option to color 1% or 5% or any...)
[similarRGB, sim_cmap, sim_labels] = similarRGBs(maskRGB,similarSnow,similarCloud,labels,cmap);
numLabels = length(sim_labels);

figure
imshow(similarRGB,sim_cmap)
h2 = colorbar;
h2.Ticks = 0.5:1:numLabels-0.5;
h2.TickLabels = sim_labels';





%% RUN ALL SCENES - USGS
for k = 1:2
    
    %similarity calculations
    if k == 1
        % 5%
        similarityFolder = 'C:\raid\scratch\timbo\snowCloud2017\similarityCalcs\fivePercent\USGS';
     per = ' 5 percent sim thresh';
    else
        % 1%
        similarityFolder = 'C:\raid\scratch\timbo\snowCloud2017\similarityCalcs\onePercent\USGS';
   
     per = ' 1 percent sim thresh';
    end




%truth masks
truthMaskFolder = 'C:\raid\data\tstillinger\snowCloud\truthMasks_matlab';

%landsat data
tiffFolder = 'C:\raid\data\tstillinger\snowCloud\snowCloudTruthMasks_Raw\ManualMasks_USGSanalyist\l8cloudmasks\clouds_snow_ice_masks';


%strings for how data is labeled
truthEnd = '_USGS_truthMasks.mat';
tiffEnd = 'data.tif';

tiffDir = dir(tiffFolder);
D = dir(similarityFolder);

for i = 1:length(D)
    if  strcmp(D(i).name(1),'L')
        
        
        %colored image of truth data with similarites
        similarityData = fullfile(similarityFolder,D(i).name);
        sceneInfo = infoFromLandsatFilename( similarityData );
        sceneID = char(sceneInfo.sceneID);
        truthData = fullfile(truthMaskFolder,[sceneID truthEnd]);
        load(similarityData)
        load(truthData)
        truthMaskSize=size(truthMask.snow);
        
        [similarSnow, similarCloud] = similarMasks(snowIdx,cloudIdx,cloudPxlSpec,snowPxlSpec,truthMaskSize);
        [ maskRGB, cmap, labels ] = maskRGBs( truthMask );
        [similarRGB, sim_cmap, sim_labels] = similarRGBs(maskRGB,similarSnow,similarCloud,labels,cmap);
        numLabels = length(sim_labels);
        
        
        %fRGB image of landsat scene
        for j = 1:length(tiffDir)
            if startsWith(tiffDir(j).name,sceneID) && endsWith(tiffDir(j).name,tiffEnd)
                [ Rtoa ] = USGSmaskTOA( fullfile(tiffFolder,tiffDir(j).name) );
                [ fRGB ] = landsat8fRGB( Rtoa );
                break
            end
        end
        
        
        %plot fRGB image next to the truth mask and similarity data
        figure
        subplot(1,2,1)
        imshow(fRGB)
        subplot(1,2,2)
        imshow(similarRGB,sim_cmap)
        h2 = colorbar;
        h2.Ticks = 0.5:1:numLabels-0.5;
        h2.TickLabels = sim_labels';
        
        suptitle([sceneID per])
        
    end
end


end



%% RUN ALL SCENES - Natalies Scenes

for k = 1:2
    if k == 1
        % 5%
        similarityFolder = 'C:\raid\scratch\timbo\snowCloud2017\similarityCalcs\fivePercent\NAT';
        per = ' 5 percent sim thresh';
    else
        % 1%
        similarityFolder = 'C:\raid\scratch\timbo\snowCloud2017\similarityCalcs\onePercent\NAT';
        per = ' 1 percent sim thresh';
    end
    
    
    %truth masks
    truthFolder = 'C:\raid\data\tstillinger\snowCloud\truthMasks_matlab\nats';
    
    %landsat data
    tiffFolder = 'C:\raid\scratch\timbo\LS8HandPaint\finalScenes';
    
    %strings for how data is labeled
    truthEnd = '_truthMasks.mat';
    tiffEnd = 'data.tif';
    
    tiffDir = dir(tiffFolder);
    D = dir(similarityFolder);
    
    for i = 1:length(D)
        if  strcmp(D(i).name(1),'L')
            
            
            %colored image of truth data with similarites
            similarityData = fullfile(similarityFolder,D(i).name);
            sceneInfo = infoFromLandsatFilename( similarityData );
            sceneID = char(sceneInfo.sceneID);
            truthData = fullfile(truthFolder,[sceneID truthEnd]);
            load(similarityData)
            load(truthData)
            truthMaskSize=size(truthMask.snow);
            
            [similarSnow, similarCloud] = similarMasks(snowIdx,cloudIdx,cloudPxlSpec,snowPxlSpec,truthMaskSize);
            [ maskRGB, cmap, labels ] = maskRGBs( truthMask );
            [similarRGB, sim_cmap, sim_labels] = similarRGBs(maskRGB,similarSnow,similarCloud,labels,cmap);
            numLabels = length(sim_labels);
            
            filepath = fullfile(tiffFolder,sceneID);
            X = GetLandsat8( filepath, 'allBands' );
            Rtoa = X.OLI_30m;
            Rtoa(:,:,9:10)=[];
            
            %rgbs from masks
            [ fRGB ] = landsat8fRGB(Rtoa);
            
            %Crop image to truth masks pixel area with 10pix buffer
            imageArea = truthMask.snow | truthMask.cloud;
            [row,col] = find(imageArea);
            ymin = min(row)-10;
            xmin = min(col)-10;
            height = max(row)+10 - ymin;
            width = max(col)+10 - xmin;
            cropRect = [xmin ymin width height];
            
            fRGB = imcrop(fRGB,cropRect);
            similarRGB = imcrop(similarRGB,cropRect);
            
            
            figure
            subplot(1,2,1)
            imshow(fRGB)
            subplot(1,2,2)
            imshow(similarRGB,sim_cmap)
            h2 = colorbar;
            h2.Ticks = 0.5:1:numLabels-0.5;
            h2.TickLabels = sim_labels';
            suptitle([sceneID per])
            
        end
    end
    
end


%%
%colormap(custom_colormap);

%colorbar





RsHomeDir = 'C:\raid\scratch\timbo\LS8HandPaint\noTruthMask\finalAGU';
outputFolder='C:\raid\scratch\timbo\LS8HandPaint\output\figures4poster\browseImgs';
landsatHomeDir = 'C:\raid\scratch\timbo\LS8HandPaint\noTruthMask\finalAGU';


% landsatFolder = 'C:\raid\scratch\timbo\LS8HandPaint\noTruthMask\LC80420342013140LGN01';

D = dir(landsatHomeDir);
for i = 1:length(D)
    if  strcmp(D(i).name(1),'L')
        similarityData = fullfile(landsatHomeDir,D(i).name);
        sceneInfo = infoFromLandsatFilename( similarityData );
        sceneID = char(sceneInfo.sceneID);
        rgbFig( similarityData,RsHomeDir,sceneID, outputFolder );
    end
end
