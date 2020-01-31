function [ Tbfix, Tb ] = testBT( )
%test function for adjustBrightnessTemp

band = [20 22 23 31 32];
solarZ = 33.33;
solarR = 579;
snowEmiss = 1-SnowCloudSensorReflectance(cosd(solarZ),250,'um','modis',band,...
    'ignoreSolar',true);

% snow fraction
mod2 = 'c:\users\dozie\Box Sync\Data\MODIS';
f = dir(fullfile(mod2,'*A2016*Endm*.h5'));
for k=1:length(f)
    modscagFile{k} = fullfile(mod2,f(k).name);
end

Tb = cell(5,1);
for k=1:5
    switch band(k)
        case 20
            Tb{k} = [269.6705 271.4360 268.3028;...
                273.9358 273.2402 270.4920;...
                271.7173 273.2402 272.0044];
        case 22
            Tb{k} = [268.0185 270.2665 266.6831;...
                272.2352 271.9515 269.2814;...
                270.7216 271.9515 270.7666];
        case 23
            Tb{k} = [266.7545 268.8743 266.2625;...
                271.1647 270.7559 267.9746;...
                269.3859 270.7559 269.5566];
        case 31
            Tb{k} = [267.4165 268.6959 265.5544;...
                270.9360 270.5585 268.4221;...
                269.4008 270.5585 269.6479];
        case 32
            Tb{k} = [267.3771 268.5476 264.8500;...
                271.1403 270.2985 268.0072;...
                269.2095 270.2985 269.5505];
        otherwise
            error('band %d not recognized',band(k))
    end
end

emiss = snowEmiss;
Tbfix = cell(5,1);
for k=1:5
    Tbfix{k} = adjustBrightnessTemp('modis',band(k),Tb{k},emiss(k),solarR);
end

end