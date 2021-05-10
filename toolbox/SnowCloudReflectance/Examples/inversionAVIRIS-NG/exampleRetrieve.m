function X = exampleRetrieve(refl,wave,R0,Z,k,m)
% X = exampleRetrieve(refl,wave,R0,Z,k)
%Example of retrievals
plotDim = [1 2];
allMethod = {'lsqnonlin','specAngle'};

% get the grain size without the dust or soot
t = wave>=900;
method = allMethod{m};

if strncmpi(method,'lsq',3)
    unk = {'radius','fsca','corrFactor'};
else
    unk = {'radius','fsca'};
end
[Ag,Bg,Cg] = invertSnowCloudSpectralRefl(refl(t,k),unk,'snow',...
    'cosZ',cosd(Z(k)),'wavelength',wave(t),'waveu','nm','R0',R0,...
    'solutionMethod',method,'waveR0',wave);
if isfield(Ag,'corrFactor')
    Cg.corrFactor = Ag.corrFactor;
end

% get the soot/dust given the grain size
[Ad,Bd,Cd] = invertSnowCloudSpectralRefl(refl(:,k),{'dust','dustr','fsca'},'snow',Cg,...
    'wavelength',wave,'R0',R0,...
    'solutionMethod',method,'waveR0',wave);

[As,Bs,Cs] = invertSnowCloudSpectralRefl(refl(:,k),{'soot','sootr','fsca'},'snow',Cg,...
    'wavelength',wave,'R0',R0,...
    'solutionMethod',method,'waveR0',wave);

% broadband albedo, for mid-latitude winter, 5 km elev, this Zenith
SP = defaultSMARTSinput('mlw','cosZ',cosd(Z(k)),'altit',4.8609);
oS = SMARTS295Main(getSMARTShome,SP);
wT = oS.spectralTbl;
radTblg = SnowCloudIntgRefl(true,wT.waveL,'nm',[wT.HorzDirect wT.HorzDiffuse],...
    'snow','bandPass',[280 4000],'waveu','nm','radius',Cg.iceRadius,...
    'cosZ',cosd(Z(k)));
radTbld = SnowCloudIntgRefl(true,wT.waveL,'nm',[wT.HorzDirect wT.HorzDiffuse],...
    'snow','bandPass',[280 4000],'waveu','nm','radius',Cg.iceRadius,...
    'cosZ',cosd(Z(k)),'dust',Ad.dust,'dustRadius',Ad.dustRadius);
radTbls = SnowCloudIntgRefl(true,wT.waveL,'nm',[wT.HorzDirect wT.HorzDiffuse],...
    'snow','bandPass',[280 4000],'waveu','nm','radius',Cg.iceRadius,...
    'cosZ',cosd(Z(k)),'soot',As.soot,'sootRadius',As.sootRadius);

X.Ag = Ag;
X.Bg = Bg;
X.Cg = Cg;
X.albedo_g = radTblg.reflectance;
X.Ad = Ad;
X.Bd = Bd;
X.Cd = Cd;
X.albedo_d = radTbld.reflectance;
X.delta_D = X.albedo_g-X.albedo_d;
X.As = As;
X.Bs = Bs;
X.Cs = Cs;
X.albedo_s = radTbls.reflectance;
X.delta_S = X.albedo_g-X.albedo_s;

%plots
subplot(plotDim(1),plotDim(2),m)
% just grain size
Rg = SnowCloudSpectralRefl(Cg,'wavelength',wave,'frac',[1 0]);
plot(wave,refl(:,k),'lineWidth',2,'LineStyle','-.','Color','k')
thisLegend = {'raw pixel'};
hold on
if isfield(Ag,'corrFactor') && ~isempty(Ag.corrFactor)
    thisLegend = cat(2,thisLegend,{'solved pixel'});
    if Bd.resnorm<Bs.resnorm
        % smooth dusty residual
        Frm = fit(wave,Bd.residual,'smoothingspline','SmoothingParam',0.009600538927466323);
        plot(wave,refl(:,k).*cosd(Z(k))/Ag.corrFactor+Frm(wave),'Color','r','LineWidth',1.5)
    else
        plot(wave,refl(:,k).*cosd(Z(k))/Ag.corrFactor+Bs.residual,'Color','r','LineWidth',1.5)
    end
elseif isfield(Bd,'spectralAngle')
    if Bd.spectralAngle<Bs.spectralAngle
        Rsolve = SnowCloudSpectralRefl(Cd);
    else
        Rsolve = SnowCloudSpectralRefl(Cs);
    end
    plot(wave,Rsolve.refl,'Color','r','LineWidth',1.5);
    thisLegend = cat(2,thisLegend,{'solved pixel'});
end

Rd = SnowCloudSpectralRefl(Cd,'frac',[1 0]);
Rs = SnowCloudSpectralRefl(Cs,'frac',[1 0]);
if isfield(Bd,'resnorm')
    X.Bd.spectralAngle = acosd(dot(Rd.refl,refl(:,k))/...
        (norm(Rd.refl)*norm(refl(:,k))));
    X.Bs.spectralAngle = acosd(dot(Rs.refl,refl(:,k))/...
        (norm(Rs.refl)*norm(refl(:,k))));
    if Bd.resnorm<Bs.resnorm
        % smooth the dusty refl
        Fsm = fit(wave,Rd.refl,'smoothingspline','SmoothingParam',0.0058451011825547);
        plot(wave,Fsm(wave),'Color',[0.93 0.69 0.13],'LineWidth',2)
        thisLegend = cat(2,thisLegend,{'dusty snow'});
    else
        plot(wave, Rs.refl,'Color',[0.64 0.078 0.18],'LineWidth',2)
        thisLegend = cat(2,thisLegend,{'sooty snow'});
    end
elseif isfield(Bd,'spectralAngle')
    if Bd.spectralAngle<Bs.spectralAngle
        plot(wave,Rd.refl,'Color',[0.93 0.69 0.13],'LineWidth',2)
        thisLegend = cat(2,thisLegend,{'dusty snow'});
    else
        plot(wave, Rs.refl,'Color',[0.64 0.078 0.18],'LineWidth',2)
        thisLegend = cat(2,thisLegend,{'sooty snow'});
    end
end
plot(wave,Rg.refl,'Color',[0 0.45 0.74],'LineWidth',2)
plot(wave,R0,'Color',[0.47 0.67 0.188],'lineWidth',2)
thisLegend = cat(2,thisLegend,{'clean snow','background'});
xlim([300 2600])
xlabel('wavelength, nm')
ylabel('spectral albedo');
if m==1
    legend(thisLegend)
end
title(method)
end