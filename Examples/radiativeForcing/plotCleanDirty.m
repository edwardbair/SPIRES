function [X] = plotCleanDirty(contam,sensor,sizeUnit,waveUnit)
% [X] = plotCleanDirty(contam,sensor,sizeUnit,waveUnit)
%make a plot of clean vs dirty snow either for spectrometer or
%multispectral sensor
%return values in a structure, depending on inputs

% set up the inputs
S = SnowCloudLimits;
r = linspace(sqrt(50),sqrt(1500),11).^2;
switch contam
    case 'dust'
        contamR = rand(1)*(S.dustRadius(2)-S.dustRadius(1))+S.dustRadius(1);
        conc = linspace(1e-4,S.dust(2),4);
    case 'soot'
        contamR = rand(1)*(S.sootRadius(2)-S.sootRadius(1))+S.sootRadius(1);
        conc = linspace(100e-9,S.soot(2),4);
    otherwise
        error('contam ''%s'' not recognized',contam)
end

% clean snow
for k=1:length(r)
    if k==1
        [R,Pc] = SnowCloudSpectralRefl('snow','cosZ',cosd(50),'sensor',sensor,...
            'waveUnit',waveUnit,'radius',r(k));
        refl = zeros(length(Pc.wavelength),length(r));
    else
        [R,Pc] = SnowCloudSpectralRefl(Pc,'radius',r(k));
    end
    refl(:,k) = R.refl;
end
cleanR = refl;
% dirty snow
% first with minimum radius, then with max
for k=1:length(conc)
    rname = [contam 'Radius'];
    if k==1
        [R,Pd] = SnowCloudSpectralRefl('snow','cosZ',cosd(50),'sensor',sensor,...
            'waveUnit',waveUnit,'radius',min(r),...
            contam,conc(k),rname,contamR);
        refl = zeros(length(Pd.wavelength),length(conc)*2);
    else
        [R,Pd] = SnowCloudSpectralRefl(Pd,contam,conc(k));
    end
    refl(:,k) = R.refl;
end
for k=1:length(conc)
    rname = [contam 'Radius'];
    if k==1
        [R,Pd] = SnowCloudSpectralRefl('snow','cosZ',cosd(50),'sensor',sensor,...
            'waveUnit',waveUnit,'radius',max(r),...
            contam,conc(k),rname,contamR);
    else
        [R,Pd] = SnowCloudSpectralRefl(Pd,contam,conc(k));
    end
    refl(:,k+length(conc)) = R.refl;
end
dirtyR = refl;

% absorption coefficient of ice
[N,w] = RefractiveIndex([],'ice',waveUnit);
k1 = find(w<min(Pd.wavelength),1,'last');
k2 = find(w>max(Pd.wavelength),1,'first');
w = w(k1:k2);
N = N(k1:k2);

% plot
ax = gca;
hdirty=plot(Pd.wavelength,dirtyR,'r','LineWidth',1.5);
hold on;
hclean=plot(Pc.wavelength,cleanR,'b','LineWidth',1);
ylabel('spectral albedo')
u = '{\mu}m';
if contains(waveUnit,'um')
    xlabel(['wavelength, ' u])
else
    xlabel(['wavelength, ' waveUnit])
end
xlim([convertLengthUnits(300,'nm',waveUnit) convertLengthUnits(2600,'nm',waveUnit)])
yyaxis right
ax.YScale = 'log';
habs=plot(w,imag(N),'k','LineWidth',2);
hglines = [hclean(1) hdirty(1) habs];
ylabel('absorption coefficient')

% cleanLabel = ['clean snow, r = ' num2str(round(r(1))) ' to ' num2str(round(r(end))) ' ' u];
if contains(sizeUnit,'um')
    cleanLabel = sprintf('clean snow, r = %d to %d %s',...
        round(r(1)), round(r(end)),u);
else
    r = convertUnits(r,'mum',sizeUnit);
    cleanLabel = sprintf('clean snow, r = %g to %g %s',...
        r(1), r(end),sizeUnit);
end
if contains(contam,'dust')
    cLabel = sprintf('dusty snow, %d to %d ppmw at min & max r',...
        round(conc(1)*1e6),round(conc(end)*1e6));
elseif contains(contam,'soot')
    cLabel = sprintf('sooty snow, %d to %d ppbw at min & max r',...
        round(conc(1)*1e9),round(conc(end)*1e9));
else
    cLabel = '';
end
absLabel = 'ice absorption coefficient';
legend(hglines,cleanLabel,cLabel,absLabel,'Location','Best')


%output
X.clearRefl = cleanR;
X.dirtyRefl = dirtyR;
X.Pc = Pc;
X.Pd = Pd;
X.radius = r;
X.concentration = conc;
end