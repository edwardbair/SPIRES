function [R1,R2] = DirectDiffuseIllustration()
%graph of direct vs diffuse illumination to show problems of identifying
%snow shadowed by topography
%

%parameters
cosZ = cosd(60);
elevation = 3000;
atmosphere = 'mlw';
grainRadius = 250;
dustConc = 50e-6;

P = defaultSMARTSinput(atmosphere,'altit',elevation,'cosZ',cosZ);
S = SMARTSMain(P);
ST = S.spectralTbl;
k1 = find(ST.waveL<400,1,'last');
k2 = find(ST.waveL>2500,1,'first');
ST = ST(k1:k2,:);
P1 = fwdPrescription('snow','radius',grainRadius,'wavelength',ST.waveL,...
    'waveU','nm','LAP','dust','LAPfraction',dustConc,'cosZ',cosZ,'lookup',false);
[~,R1] = SPIReS_fwd(P1);
P2 = fwdPrescription('snow','radius',grainRadius,'wavelength',ST.waveL,...
    'waveU','nm','LAP','dust','LAPfraction',dustConc,'cosZ',[],'lookup',false);
[~,R2] = SPIReS_fwd(P2);

%graph
Ftot = fit(ST.waveL,(ST.HorzDirect.*R1.reflectance+ST.HorzDiffuse.*R2.reflectance)./...
    (ST.HorzDirect + ST.HorzDiffuse),'smoothingspline');
Fdif = fit(ST.waveL,(ST.HorzDiffuse.*R2.reflectance)./...
    (ST.HorzDirect + ST.HorzDiffuse),'smoothingspline');
% add MODIS bands
hold on;
MT = SensorTable('modis','nm');
grayColor = [.75 .75 .75];
for b=1:7
    w = MT.UpperWavelength(b)-MT.LowerWavelength(b);
    ht = Ftot(MT.CentralWavelength(b));
    if b==1
        h1 = bar(MT.CentralWavelength(b),ht,w,'EdgeColor',grayColor,'FaceColor',grayColor);
    else
        bar(MT.CentralWavelength(b),ht,w,'EdgeColor',grayColor,'FaceColor',grayColor)
    end
end
% plot spectra
h2 = plot(ST.waveL,[Ftot(ST.waveL) Fdif(ST.waveL)],'LineWidth',2);
xlim([300 2600])
xlabel('wavelength, nm')

legend([h1; h2],{'MODIS bands','snow reflectance','reflectance of snow in topographic shade'})
end
