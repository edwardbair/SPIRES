function [] = inversionGraphsWRR(truthFilename,LS8C1_folder,whichData,pixelT)
%inversion snow cloud matches for WRR paper

%LOAD DATA
fileInfo = infoFromLandsatFilename(truthFilename);
sceneID = char(fileInfo.sceneID);

switch whichData
    case 'usgs'
        
        geotiffFilename = truthFilename;
        geotiffFilename(end-7:end)='data.tif';
        truthSceneInfo = geotiffinfo(geotiffFilename);
        x = truthSceneInfo.CornerCoords.X;
        y = truthSceneInfo.CornerCoords.Y;
        
    case 'nats'
        
        load(truthFilename)
        %crop the truth mask to just the region of the image with the truth pixels in it
        imageArea = truthMask.snow | truthMask.cloud |truthMask.empty;
        [row,col] = find(imageArea);
        
        rowmin = min(row);
        colmin = min(col);
        rowmax = max(row);
        colmax = max(col);
        rowBox = [rowmin rowmax rowmax rowmin];
        colBox = [colmin colmin colmax colmax];
        
        adjsutRow = [-.5 .5 .5 -.5];
        adjustCol = [-.5 -.5 .5 .5];
        rowBox = rowBox+adjsutRow;
        colBox = colBox + adjustCol;
        
        R= truthMask.sceneInfo.RefMatrix;
        [x,y] = pix2map(R,rowBox,colBox);
        
end

% LOAD LANDSAT 8 COLLECTION 1 DATA
X = GetLandsat8_Collection1( LS8C1_folder, 'allbands' );
[~,~,fullSceneInfo] = GetLandsat8_Collection1( LS8C1_folder, 'BQA' );

LS8_30m = X.OLI_30m;

% FIND TRUTH MASK IN FULL LANDSAT SCENE
[Rtoa,~] = extractSubscene(fullSceneInfo.RefMatrix,x,y,LS8_30m,whichData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Visualize Results %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fRGBbox = landsat8fRGB(Rtoa, [6,5,4]);

pxlIdx = false(size(Rtoa));
pxlIdx(pixelT.pixelIDX) = true;

[r,c]=find(pxlIdx);
testPxl = false(size(101,101));
testPxl(50,50)=true;
rmin = r-50; rmax = r+50;
cmin = c-50; cmax = c+50;
subFRGB = fRGBbox(rmin:rmax,cmin:cmax,:);

%OLI spectra graph
sensor='landsatOLI';
measured = pixelT.measuredReflectance;
cloudModel =measured(1:7) + pixelT.cloudResidual;
snowModel =measured(1:7) + pixelT.snowResidual;
bandWvs = SensorTable(sensor);
bandWvs = bandWvs.CentralWavelength(1:7);

% hyperspectral model of the snow and cloud inversion
hyperSpec = SensorTable('AVIRIS-NG');
lambda = hyperSpec.CentralWavelength;
%other fractions
soilType = 'darkLoam';
soilR = SoilReflectance( lambda, soilType);
vegType = 'conifer';
vegR = VegetationReflectance( lambda, vegType );
ZA = pixelT.solarZ;
snowRadius = pixelT.Sradius;
cloudRadius = pixelT.Cradius;
waveUnit = 'um';

switch char(pixelT.Background)
    case 'Dark Loam + shade'       
        R0 = table(lambda,soilR,'VariableNames',{'wavelength','reflectance'});
        R0.reflectance = [R0.reflectance zeros(height(R0),1)];
        otherEndMem = pixelT.otherEndMem;
        
    case 'Conifer'
        R0 = table(lambda,vegR,'VariableNames',{'wavelength','reflectance'});
        otherEndMem = pixelT.otherEndMem(1);
    
    case 'Dark Loam'
        R0 = table(lambda,soilR,'VariableNames',{'wavelength','reflectance'});
        otherEndMem = pixelT.otherEndMem(1);
        
end

R0.Properties.VariableUnits = {'nm',''};

% % % pixelAlt = 3; 
% % % APS = defaultSMARTSinput('mlw','cosZ',cosd(ZA),'altit',pixelAlt);
% % % SS = SMARTS295Main(getSMARTShome,APS);
% % % solarTbl = table(SS.spectralTbl.waveL,...
% % %     [SS.spectralTbl.HorzDirect SS.spectralTbl.HorzDiffuse],...
% % %     'VariableNames',{'wavelength','irradiance'});
% % % solarTbl.Properties.VariableUnits = {'nm','W/m2/nm'};

%Jeff's inversions output but with spectrometer
snowR = SnowCloudSpectralRefl('snow',...
    'R0',R0,'radius',snowRadius,'dust',pixelT.dust,'cosZ',cosd(ZA),...
    'fractionalCoverage',[pixelT.fSCA otherEndMem],'waveUnit',waveUnit,'wavelength',R0.wavelength);

cloudR = SnowCloudSpectralRefl('mixedcloud',...
    'R0',R0,'radius',cloudRadius,'cosZ',cosd(ZA),'waterEquivalent',...
    pixelT.waterEq,'waveUnit',waveUnit,'wavelength',R0.wavelength);

% % % % SNOWCLOUDINTGREFL
% % % snowT = SnowCloudIntgRefl(solarTbl,'snow',...
% % %     'R0',R0,'radius',snowRadius,'dust',pixelT.dust,'cosZ',cosd(ZA),'fractionalCoverage',[pixelT.fSCA otherEndMem],'sensor','AVIRIS-NG');
% % % 
% % % cloudT = SnowCloudIntgRefl(solarTbl,'mixedCloud',...
% % %     'R0',R0,'radius',cloudRadius,'cosZ',cosd(ZA),'sensor','AVIRIS-NG','waterEquivalent',pixelT.waterEq);
% % % 

% plot
figure
subplot(1,2,1)
imshow(subFRGB)
hold on
SE = strel('square',10);
visboundaries(imdilate(testPxl,SE),'Color','r');
xlabel(sceneID)
set(gca,'FontSize',20)

% % subplot(1,3,2)
% % scatter(bandWvs.*1000,measured(1:7),100,'+','LineWidth',2)
% % hold on
% % scatter(bandWvs.*1000,cloudModel,100,'*','LineWidth',2)
% % scatter(bandWvs.*1000,snowModel,100,'x','LineWidth',2)
% % ylim([0 1])
% % xlim([350 2500])
% % legend('measured','cloud model','snow model')
% % axis('square')
% % set(gca,'FontSize',20)
% % ylabel('spectral albedo')
% % xlabel('wavelength (nm)')

subplot(1,2,2)
olibandpass = SensorTable('landsatoli');
for k = 1:7
h(k) = area([olibandpass.LowerWavelength(k) olibandpass.UpperWavelength(k)],[measured(k) measured(k)]);
h(k).FaceColor = [0.8 0.8 0.8];
h(k).EdgeAlpha=0;
if k ==1 
    hold on
end

hsc = plot(R0.wavelength,cloudR.refl,'r','LineWidth',1.5);
hss = plot(R0.wavelength,snowR.refl,'b','LineWidth',1.5);

% hss = plot(mean(snowT.bandPass,2),snowR.refl,'b','LineWidth',1.5);
%ms = scatter(bandWvs.*1000,measured(1:7),100,'k','+','LineWidth',2);

msc = scatter(bandWvs,cloudModel,100,'r','*','LineWidth',2);
mss = scatter(bandWvs,snowModel,100,'b','x','LineWidth',2);

ylabel('spectral albedo')
xlabel('wavelength (um)')
ylim([0 1])
xlim([.350 2.5])
axis('square')
legend([h(1) msc mss hsc hss ],'measured OLI bandpasses','cloud model OLI','snow model OLI','cloud model spectrometer','snow model spectrometer')
set(gca,'FontSize',20)
end


