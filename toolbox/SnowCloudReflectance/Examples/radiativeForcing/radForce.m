function [outS] = radForce(snowSize,sootConc)
% [outS] = radForce(snowSize,sootConc)
%graph to illustrate radiative forcing
%

% input spectrum from SMARTS
cosZ = cosd(60);
altit = 3; % km
P = defaultSMARTSinput('mlw','cosZ',cosZ,'altit',altit);
S = SMARTS295Main(getSMARTShome,P);
sTbl = S.spectralTbl;
outS.waveL = sTbl.waveL;
outS.irradiance = [sTbl.HorzDirect sTbl.HorzDiffuse];

% snow albedo
Rc = SnowCloudSpectralRefl('snow','radius',snowSize,'sizeu','um',...
    'wavelength',outS.waveL,'waveu','nm','cosZ',cosZ);
Rcd = SnowCloudSpectralRefl('snow','radius',snowSize,'sizeu','um',...
    'wavelength',outS.waveL,'waveu','nm','cosZ',[]);
Rd = SnowCloudSpectralRefl('snow','radius',snowSize,'sizeu','um',...
    'wavelength',outS.waveL,'waveu','nm','cosZ',cosZ,'soot',sootConc);
Rdd = SnowCloudSpectralRefl('snow','radius',snowSize,'sizeu','um',...
    'wavelength',outS.waveL,'waveu','nm','cosZ',[],'soot',sootConc);
outS.reflC = [Rc.refl Rcd.refl];
outS.reflD = [Rd.refl Rdd.refl];

% reflected radiation
upRadC = outS.irradiance(:,1).*Rc.refl + outS.irradiance(:,2).*Rcd.refl;
upRadD = outS.irradiance(:,1).*Rd.refl + outS.irradiance(:,2).*Rdd.refl;
outS.upRadC = upRadC;
outS.upRadD = upRadD;

% net radiation
netRadC = sum(outS.irradiance,2)-outS.upRadC;
netRadD = sum(outS.irradiance,2)-outS.upRadD;
outS.netRadC = netRadC;
outS.netRadD = netRadD;

end