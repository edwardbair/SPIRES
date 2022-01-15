function [T,F,G] = tblSizeIndices()
%size index for grains
%

wv = [985 1035 1090];
rad = linspace(sqrt(50),sqrt(1000),50).^2;
cosZ = cosd(60);
[r,w] = ndgrid(rad,wv);
% clean snow
prescription = fwdPrescription('snow','radius',r(:),'wavelength',w(:),'waveunit','nm',...
    'cosZ',cosZ);
[~,T] = SPIReS_fwd(prescription);
t1 = T.wavelength==wv(1);
t2 = T.wavelength==wv(2);
t3 = T.wavelength==wv(3);
y = sqrt(T.radius(t2));

%fits
[Fboth,Gboth] = fit((mean([T.reflectance(t1) T.reflectance(t3)],2)-...
    T.reflectance(t2))./mean([T.reflectance(t1) T.reflectance(t3)],2),...
    y,'power2','robust','bisquare');
[F985,G985] = fit((T.reflectance(t1)-T.reflectance(t2))./T.reflectance(t1),...
    y,'power2','robust','bisquare');
[F1090,G1090] = fit((T.reflectance(t3)-T.reflectance(t2))./T.reflectance(t3),...
    y,'power2','robust','bisquare');

F{1} = F985;
F{2} = F1090;
F{3} = Fboth;

G{1} = G985;
G{2} = G1090;
G{3} = Gboth;

%plots
for k=1:3
    switch k
        case 1
            x = (T.reflectance(t1)-T.reflectance(t2))./T.reflectance(t1);
            hold on;
        case 2
            x = (T.reflectance(t3)-T.reflectance(t2))./T.reflectance(t3);
        case 3
            x = (mean([T.reflectance(t1) T.reflectance(t3)],2)-...
                T.reflectance(t2))./mean([T.reflectance(t1) T.reflectance(t3)],2);
    end
    scatter(x,T.radius(t2),'.')
    hold on;
    plot(x,(F{k}(x)).^2,'linewidth',1)
end
ylabel('snow radius, {\mu}m')
xlabel('index (see text)')
end