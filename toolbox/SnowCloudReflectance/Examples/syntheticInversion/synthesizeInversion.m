function X=synthesizeInversion(invMeth)
% X=synthesizeInversion(invMeth)
%% Title
% Synthesize spectra, add noise, and try to retrieve inputs

%% input values
wavelength = linspace(380,2500,425)';
noise = 0.01;
bias = 0;
nValues = 500;
S = SnowCloudLimits;
Nradius = 50;
botratio = 5/3;
topratio = 12/20
Uradius = linspace(sqrt(S.snowRadius(1)*botratio),sqrt(S.snowRadius(2)*topratio),Nradius).^2;
NfSCA = 50;
UfSCA = linspace(0.01,1,NfSCA);
Uelev = 2:.5:6;
UsolarZ = 30:10:60;
Ndust = 10;
Nsoot = 10;
UdustRc = logspace(log10(S.dustRadius(1)*botratio),log10(S.dustRadius(2)*topratio),Ndust);
UdustC = linspace(S.dust(1),S.dust(2),Ndust);
UsootRc = logspace(log10(S.sootRadius(1)*botratio),log10(S.sootRadius(2)*topratio),Nsoot);
UsootC = linspace(S.soot(1),S.soot(2),Nsoot);
%% all possible combinations of inputs
[ssa,fSCA,elev,solarZ,dustRc,dustC] = ndgrid(radius2SSA(Uradius,S.unitsSize),...
    UfSCA,Uelev,UsolarZ,UdustRc,UdustC);
[~,~,~,~,sootRc,sootC] = ndgrid(radius2SSA(Uradius,S.unitsSize),...
    UfSCA,Uelev,UsolarZ,UsootRc,UsootC);
ssa = ssa(:);
fSCA = fSCA(:);
elev = elev(:);
solarZ = solarZ(:);
dustRc = dustRc(:);
dustC = dustC(:);
sootRc = sootRc(:);
sootC = sootC(:);
inputTbl = table(ssa,fSCA,elev,solarZ,dustRc,dustC,sootRc,sootC);
nTotal = height(inputTbl);
% randomly select nValues
rng('shuffle')
kSelect = randperm(nTotal,nValues);
inputTbl = inputTbl(kSelect',:);

%% for each row in inputTbl, generate dust and soot spectra, add noise and
% bias, and then try to retrieve input properties
% input spectrum for each elevation needed
atmos = 'mlw';
R0 = SoilReflectance(wavelength,'DarkLoam','nm');
T = inputTbl; % (just to shorten the commands)
tradwave = wavelength>=900;
for k=1:height(T)
    [Rd,P] = SnowCloudSpectralRefl('snow','ssa',T.ssa(k),'fSCA',T.fSCA(k),...
        'cosZ',cosd(T.solarZ(k)),'dustradius',T.dustRc(k),'dust',T.dustC(k),...
        'wavelength',wavelength,'waveu','nm','R0',R0,'waveR0',wavelength);
    Rs = SnowCloudSpectralRefl(P,'sootradius',T.sootRc(k),'soot',T.sootC(k));
    % add noise
    n1 = randn(length(wavelength),1)*noise+bias;
    n2 = randn(length(wavelength),1)*noise+bias;
    rd = Rd.refl.*(1+n1);
    rs = Rs.refl.*(1+n2);
    % retrieve ssa first
    [L,stat] = invertSnowCloudSpectralRefl(rd(tradwave),{'ssa','fSCA'},'snow',...
        'wavelength',wavelength(tradwave),'R0',R0(tradwave),'waveu','nm',...
        'solutionMethod',invMeth,'waveR0',wavelength(tradwave),'cosZ',cosd(T.solarZ(k)));
    if contains(invMeth,'lsq','IgnoreCase',true)
        thisTbl = table(L.ssa,stat.resnorm,'VariableNames',{'ssa','resnorm'});
    else
        thisTbl = table(L.ssa,stat.spectralAngle,'VariableNames',{'ssa','spectralAngle'});
    end
    % now retrieve dust and soot
    [LD,statD,Px] = invertSnowCloudSpectralRefl(rd,{'fSCA','dust','dustradius'},'snow',...
        'wavelength',wavelength,'R0',R0,'waveu','nm','waveR0',wavelength,...
        'solutionMethod',invMeth,'ssa',L.ssa,'cosZ',cosd(T.solarZ(k)));
    [LS,statS] = invertSnowCloudSpectralRefl(rs,{'fSCA','soot','sootradius'},'snow',...
        'wavelength',wavelength,'R0',R0,'waveu','nm','waveR0',wavelength,...
        'solutionMethod',invMeth,'ssa',L.ssa,'cosZ',cosd(T.solarZ(k)));
    thisTbl = [thisTbl table(LD.fSCA,LD.dust,LD.dustRadius,LS.fSCA,LS.soot,LS.sootRadius,...
        [statD.exitflag statS.exitflag],'VariableNames',...
        {'fSCA_dust','dustC','dustRc','fSCA_soot','sootC','sootRc','exitflags'})]; %#ok<AGROW>
    if contains(invMeth,'lsq','IgnoreCase',true)
        thisTbl = [thisTbl table(statD.resnorm,statS.resnorm,'VariableNames',...
            {'dust_resnorm','soot_resnorm'})]; %#ok<AGROW>
    else
        thisTbl = [thisTbl table(statD.spectralAngle,statS.spectralAngle,'VariableNames',...
            {'dust_specAngle','soot_specAngle'})]; %#ok<AGROW>
    end
    if k==1
        oTbl = thisTbl;
    else
        oTbl = [oTbl; thisTbl]; %#ok<AGROW>
    end
end
X.solutionMethod = Px.solutionMethod;
X.oTbl = oTbl;
X.inputTbl = T;
% compare clean albedo with dusty and sooty albedo
% all combinations of elevation and solarZ
[elev,solarZ] = ndgrid(Uelev,UsolarZ);
elev = elev(:);
solarZ = solarZ(:);
for k=1:length(elev)
    SMP = defaultSMARTSinput(atmos,'cosZ',cosd(solarZ(k)),'altit',elev(k));
    XS = SMARTS295Main(getSMARTShome,SMP);
    waveL = XS.spectralTbl.waveL;
    irrad{k} = [XS.spectralTbl.HorzDirect XS.spectralTbl.HorzDiffuse]; %#ok<AGROW>
end
for k=1:height(X.oTbl)
    m = elev==X.inputTbl.elev(k) & solarZ==X.inputTbl.solarZ(k);
    cleanA = SnowCloudIntgRefl(true,waveL,'nm',irrad{m},'snow','bandPass',...
        [280 4000],'waveu','nm','ssa',X.inputTbl.ssa(k),'cosZ',cosd(X.inputTbl.solarZ(k)));
    cleanA_solve = SnowCloudIntgRefl(true,waveL,'nm',irrad{m},'snow','bandPass',...
        [280 4000],'waveu','nm','ssa',X.oTbl.ssa(k),'cosZ',cosd(X.inputTbl.solarZ(k)));
    dustyA = SnowCloudIntgRefl(true,waveL,'nm',irrad{m},'snow','bandPass',...
        [280 4000],'waveu','nm','ssa',X.inputTbl.ssa(k),'cosZ',cosd(X.inputTbl.solarZ(k)),...
        'dust',X.inputTbl.dustC(k),'dustR',X.inputTbl.dustRc(k));
    dustyA_solve = SnowCloudIntgRefl(true,waveL,'nm',irrad{m},'snow','bandPass',...
        [280 4000],'waveu','nm','ssa',X.oTbl.ssa(k),'cosZ',cosd(X.inputTbl.solarZ(k)),...
        'dust',X.oTbl.dustC(k),'dustR',X.oTbl.dustRc(k));
    sootyA = SnowCloudIntgRefl(true,waveL,'nm',irrad{m},'snow','bandPass',...
        [280 4000],'waveu','nm','ssa',X.inputTbl.ssa(k),'cosZ',cosd(X.inputTbl.solarZ(k)),...
        'soot',X.inputTbl.sootC(k),'sootR',X.inputTbl.sootRc(k));
    sootyA_solve = SnowCloudIntgRefl(true,waveL,'nm',irrad{m},'snow','bandPass',...
        [280 4000],'waveu','nm','ssa',X.oTbl.ssa(k),'cosZ',cosd(X.inputTbl.solarZ(k)),...
        'soot',X.oTbl.sootC(k),'sootR',X.oTbl.sootRc(k));
    thisTbl = table(cleanA.reflectance, cleanA_solve.reflectance,...
        dustyA.reflectance, dustyA_solve.reflectance,...
        sootyA.reflectance, sootyA_solve.reflectance,'VariableNames',...
        {'cleanA','cleanA_solve','dustyA','dustyA_solve','sootyA','sootyA_solve'});
    if k==1
        albedoTbl = thisTbl;
    else
        albedoTbl = [albedoTbl; thisTbl]; %#ok<AGROW>
    end
end
X.albedoTbl = albedoTbl;
end
