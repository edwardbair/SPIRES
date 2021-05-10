function [] = inversionPlot(pixelT)
%inversion snow/cloud matches for WRR paper, like Timbo's code but w/o images

waveUnit = 'nm';
%OLI spectra graph
sensor='landsatOLI';
measured = pixelT.measuredReflectance;
% cloudModel = measured(1:7) + pixelT.cloudResidual;
% snowModel = measured(1:7) + pixelT.snowResidual;
% bandWvs = SensorTable(sensor,waveUnit);
% bandWvs = bandWvs.CentralWavelength(1:7);

% hyperspectral model of the snow and cloud inversion
hyperSpec = SensorTable('AVIRIS-NG',waveUnit);
lambda = hyperSpec.CentralWavelength;
%other fractions
soilType = 'darkLoam';
soilR = SoilReflectance( lambda, soilType,waveUnit);
vegType = 'conifer';
vegR = VegetationReflectance( lambda, vegType, waveUnit);
ZA = pixelT.solarZ;
snowRadius = pixelT.Sradius;
cloudRadius = pixelT.Cradius;


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

%modeled OLI values
AS = defaultSMARTSinput('mlw','cosZ',cosd(ZA));
[~,SMARTStbl] = SMARTS295Main(getSMARTShome,AS);
snowL = SnowCloudIntgRefl(SMARTStbl,'snow',...
   'R0',R0,'radius',snowRadius,'dust',pixelT.dust,'cosZ',cosd(ZA),...
    'fractionalCoverage',[pixelT.fSCA otherEndMem],'waveUnit',waveUnit,...
    'sensor',sensor,'bands',1:7);
cloudL = SnowCloudIntgRefl(SMARTStbl,'mixed',...
   'R0',R0,'radius',cloudRadius,'dust',pixelT.dust,'cosZ',cosd(ZA),...
    'waterEq',pixelT.waterEq,'waveUnit',waveUnit,...
    'sensor',sensor,'bands',1:7,'wetness',pixelT.wetness);

% % % % SNOWCLOUDINTGREFL
% % % snowT = SnowCloudIntgRefl(solarTbl,'snow',...
% % %     'R0',R0,'radius',snowRadius,'dust',pixelT.dust,'cosZ',cosd(ZA),'fractionalCoverage',[pixelT.fSCA otherEndMem],'sensor','AVIRIS-NG');
% % %
% % % cloudT = SnowCloudIntgRefl(solarTbl,'mixedCloud',...
% % %     'R0',R0,'radius',cloudRadius,'cosZ',cosd(ZA),'sensor','AVIRIS-NG','waterEquivalent',pixelT.waterEq);
% % %

% plot
% figure
% subplot(1,2,1)
% imshow(subFRGB)
% hold on
% SE = strel('square',10);
% visboundaries(imdilate(testPxl,SE),'Color','r');
% xlabel(sceneID)
% set(gca,'FontSize',20)

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

% subplot(1,2,2)
olibandpass = SensorTable(sensor,waveUnit);
for k = 1:7
    h = bar(olibandpass.CentralWavelength(k),measured(k),...
        olibandpass.UpperWavelength(k)-olibandpass.LowerWavelength(k),...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
    if k ==1
        hold on
    end
end

hsc = plot(R0.wavelength,cloudR.refl,'r','LineWidth',1.5);
hss = plot(R0.wavelength,snowR.refl,'b','LineWidth',1.5);
hsb = plot(R0.wavelength,mean(R0.reflectance,2),'k','LineWidth',1,'LineStyle','--');

% hss = plot(mean(snowT.bandPass,2),snowR.refl,'b','LineWidth',1.5);
%ms = scatter(bandWvs.*1000,measured(1:7),100,'k','+','LineWidth',2);

msc = scatter(mean(cloudL.bandPass,2),cloudL.reflectance,100,'r','*','LineWidth',2);
mss = scatter(mean(snowL.bandPass,2),snowL.reflectance,100,'b','x','LineWidth',2);

ylabel('spectral albedo')
xlabel(['wavelength (' waveUnit ')'])
ylim([0 1])
if strcmpi(waveUnit,'mum') || strcmpi(waveUnit,'um')
    xlim([.3 2.6])
elseif strcmpi(waveUnit,'nm')
    xlim([300 2600]);
end
axis('square')
legend([h msc mss hsc hss hsb ],'measured OLI reflectance','OLI cloud model',...
    'OLI snow model','spectral cloud model','spectral snow model','background')
set(gca,'FontSize',12)
ff = pixelT.file{1};
uscore = strfind(ff,'_');
title([ff(1:uscore-1) ' ' num2str(pixelT.pixelIDX)])
end