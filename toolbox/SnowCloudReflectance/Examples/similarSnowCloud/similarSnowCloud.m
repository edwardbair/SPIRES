function [cloudP,snowP] = similarSnowCloud(file,nSample,cloudType,solutionMethod)
% [cloudP,snowP] = similarSnowCloud(file,nSample,cloudType,solutionMethod)
% solve for cloud and snow properties for the same spectrum
% file - files from C:\raid\scratch\timbo\ch1_stats\dec2018\FPFNspec with
%   reflectance values in columns 1-8 (bands 1:7 9] and solar zenith angle in column 9,
% nSample - randomly select this many samples from each category
%   (FNcloudSnow, FNccloud, FPsnow)
% cloudType - 'iceCloud', 'waterCloud', or 'mixedCloud'
% solutionMethod - either 'lsqnonlin' or 'specAngle'

% output structure
cloudP = struct([]);
snowP = struct([]);

% load file, parse structure, and select nSample from each category
CF = load(file);
if nSample>=size(CF.FPFNspec.FNcloudSnow,1)
    p = randperm(size(CF.FPFNspec.FNcloudSnow,1));
else
    p = randperm(size(CF.FPFNspec.FNcloudSnow,1),nSample);
end
if isempty(p)
    FNcloudSnow = [];
else
    FNcloudSnow = CF.FPFNspec.FNcloudSnow(p',:);
    FNcloudSnow(FNcloudSnow<0) = 0;
end
p = randperm(size(CF.FPFNspec.FNcloud,1),nSample);
FNcloud = CF.FPFNspec.FNcloud(p',:);
FNcloud(FNcloud<0) = 0;
if nSample>=size(CF.FPFNspec.FPsnow,1)
    p = randperm(size(CF.FPFNspec.FPsnow,1));
else
    p = randperm(size(CF.FPFNspec.FPsnow,1),nSample);
end
FPsnow = CF.FPFNspec.FPsnow(p',:);
FPsnow(FPsnow<0) = 0;

% run through all combinations for 4 different backgrounds
for s=1:nSample
    for t=1:3
        [~,filename,~] = fileparts(file);
        cloudP(s,t).file = filename;
        snowP(s,t).file = filename;
        switch t
            case 1
                if isempty(FNcloudSnow)
                    cloudP(s,t).ErrorType = 'empty matrix';
                    thisOne = [];
                else
                    thisOne = FNcloudSnow(s,:);
                    cloudP(s,t).ErrorType = 'FNcloudSnow';
                end
            case 2
                thisOne = FNcloud(s,:);
                cloudP(s,t).ErrorType = 'FNcloud';
            case 3
                if s<=size(FPsnow,1)
                    thisOne = FPsnow(s,:);
                    cloudP(s,t).ErrorType = 'FPsnow';
                else
                    thisOne = [];
                    cloudP(s,t).ErrorType = 'end of matrix';
                end
        end
        snowP(s,t).ErrorType = cloudP(s,t).ErrorType;
        snowP(s,t).solutionMethod = solutionMethod;
        cloudP(s,t).solutionMethod = solutionMethod;
        if isempty(thisOne)
            cloudP(s,t).outputStruct = struct([]);
            snowP(s,t).outputStruct = struct([]);
            continue
        end
        cloudP(s,t).refl = thisOne(1:8);
        snowP(s,t).refl = thisOne(1:7);
        cloudP(s,t).solarZ = thisOne(end);
        snowP(s,t).solarZ = thisOne(end);
        cloudP(s,t).pathrow = [floor(thisOne(10)/100) mod(thisOne(10),100)];
        snowP(s,t).pathrow = [floor(thisOne(10)/100) mod(thisOne(10),100)];
        cloudP(s,t).dnval = thisOne(9);
        snowP(s,t).dnval = thisOne(9);
        cloudP(s,t).idx = thisOne(11);
        snowP(s,t).idx = thisOne(11);
        
        APS = defaultSMARTSinput('mlw','cosZ',cosd(thisOne(end)),'altit',3);
        [~,STbl] = SMARTS295Main(getSMARTShome,APS);
        DarkLoam = table(STbl.wavelength,...
            SoilReflectance(STbl.wavelength,'darkloam',STbl.Properties.VariableUnits{1}),...
            'VariableNames',{'wavelength','reflectance'});
        DarkLoam.Properties.VariableUnits = {STbl.Properties.VariableUnits{1},''};
        Conifer = table(STbl.wavelength,...
            VegetationReflectance(STbl.wavelength,'conifer',STbl.Properties.VariableUnits{1}),...
            'VariableNames',{'wavelength','reflectance'});
        Conifer.Properties.VariableUnits = {STbl.Properties.VariableUnits{1},''};
        for b=1:3
            switch b
                case 1
                    R0 = DarkLoam;
                    Background = 'Dark Loam';
                case 2
                    R0 = DarkLoam;
                    R0.reflectance = [R0.reflectance zeros(height(R0),1)];
                    Background = 'Dark Loam + shade';
                    %                 case 3
                    %                     R0 = [Conifer DarkLoam];
                    %                     Background = 'Conifer + Dark Loam';
                case 3
                    R0 = Conifer;
                    Background = 'Conifer';
            end
            
            % solve assuming snow
            snowStruct.Background = Background;
            try
                % solve for grain size first
                [o,ss,P] = invertSnowCloudIntgRefl(STbl,thisOne(4:7),...
                    {'radius','fSCA'},'snow','cosZ',cosd(thisOne(end)),'R0',R0,...
                    'sensor','landsatoli','bands',4:7,...
                    'solutionMethod',solutionMethod,'waveu',...
                    STbl.Properties.VariableUnits{1}); %#ok<ASGLU>
                % now use that grain size to solve for other variables
                [oStruct,stats,P] = invertSnowCloudIntgRefl(STbl,thisOne(1:7),...
                    {'fSCA','dust'},'snow','radius',o.radius,'cosZ',cosd(thisOne(end)),'R0',R0,...
                    'sensor','landsatoli','bands',1:7,...
                    'solutionMethod',solutionMethod,'waveu',...
                    STbl.Properties.VariableUnits{1});
                oStruct.radius = o.radius;
                snowStruct.oStruct = oStruct;
                snowStruct.stats = stats;
                snowStruct.P = P;
            catch
                snowStruct.oStruct = struct([]);
                snowStruct.stats = struct([]);
                snowStruct.P = struct([]);
            end
            
            % solve assuming cloud
            solveV = {'radius','waterEquivalent'};
            if contains(cloudType,'mixed')
                solveV = {'radius','waterEquivalent','wetness'};
            end
            cloudStruct.Background = Background;
            try
                [oStruct,stats,P] = invertSnowCloudIntgRefl(STbl,thisOne(1:7),...
                    solveV,cloudType,'cosZ',cosd(thisOne(end)),'R0',R0,...
                    'sensor','landsatoli','bands',1:7,...
                    'solutionMethod',solutionMethod,...
                    'waveu',STbl.Properties.VariableUnits{1});
                cloudStruct.oStruct = oStruct;
                cloudStruct.stats = stats;
                cloudStruct.P = P;
            catch
                cloudStruct.oStruct = struct([]);
                cloudStruct.stats = struct([]);
                cloudStruct.P = struct([]);
            end
            % put into output structure
            cloudP(s,t).outputStruct(b) = cloudStruct;
            snowP(s,t).outputStruct(b) = snowStruct;
        end
    end
end
end