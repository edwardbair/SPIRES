function [C,varargout]=generateBroadbandAlbedo(contamName,atmos,altit,nRadii,nAngle)
% [C [,F]]=generateBroadbandAlbedo(contamName,atmos,altit,nRadii,nAngle)
% function to compute albedo values for snow radii from 30 to 2000 um and
% cosines of solar Z from 0.05 to 1, using SMARTS input for specified atmosphere
% (e.g., 'mlw', 'sas') and altitude
% optional output, generate the griddedInterpolant lookup table

rad = linspace(sqrt(30),sqrt(2000),nRadii);
rad = round(rad.^2)';
cosZ = linspace(.05,1,round(nAngle/4))';
cosZp = linspace(.05,1,nAngle)';
if contains(contamName,'dust')
    contam = 1e-3*linspace(0,1,11);
elseif contains(contamName,'soot')
    contam = 5e-5*linspace(0,1,11);
end

% all possible SMARTS295 combinations of altit & cosZ
[aS,cS] = ndgrid(altit,cosZp);
AltCosineCombo = table(aS(:),cS(:),'VariableNames',{'altit','cosZ'});
for k=1:height(AltCosineCombo)
    atmosP = defaultSMARTSinput(atmos,'cosZ',AltCosineCombo.cosZ(k),...
        'altit',AltCosineCombo.altit(k));
    [~,T] = SMARTS295Main(getSMARTShome,atmosP);
    STbl{k} = T; %#ok<AGROW>
end
str1 = sprintf('completed SMARTS295 %d combinations of %d cosZ and %d altit',...
    length(STbl),length(cosZ),length(altit));
disp(str1) %#ok<DSPS>

% all combinations of loop variables
[loopContam,loopAltitude,loopCosine,loopCosZ] = ndgrid(contam,altit,cosZp,cosZ);
loopContam = loopContam(:);
loopAltitude = loopAltitude(:);
loopCosine = loopCosine(:);
loopCosZ = loopCosZ(:);
% matrix of reflectance values
reflMat = zeros(length(rad),length(loopCosine));

str2 = sprintf('now considering %d combinations of %d cosZ, %d cosZp, %d altit, and %d %s, with %d radii',...
    length(loopCosine),length(cosZ),length(cosZp),length(altit),length(contam),...
    contamName,length(rad));
disp('')
disp(str2) %#ok<DSPS>
parfor k=1:length(loopCosine)
    t = AltCosineCombo.cosZ==loopCosine(k) &...
        AltCosineCombo.altit==loopAltitude(k);
    kt = find(t);
    assert(length(kt)==1,...
        'SMARTSinputTbl should find just 1 match, found %d, index %d',...
        length(kt),k)
    T = STbl{kt};
    
    % multitple radii match this cosZ, altit combination
    if loopContam(k)==0
        [R,~] = SnowCloudIntgRefl(T,'snow','cosZ',loopCosine(k),'bandPass',...
            [280 4000],'radius',rad);
    else
        [R,~] = SnowCloudIntgRefl(T,'snow','cosZ',loopCosine(k),'bandPass',...
            [280 4000],'radius',rad,contamName,loopContam(k));
    end
    reflMat(:,k) = R.reflectance;
    if mod(k,1000)==0
        disp(['combo ' num2str(k)])
    end
end

% convert inputs and output into 4D arrays that griddedInterpolant can use
[radius,contamGrid,altitGrid,cosZgrid,cosZpgrid] =...
    ndgrid(rad,contam,altit,cosZ,cosZp);
refl = zeros(size(radius));
for k=1:length(loopCosine)
    t = cosZpgrid==loopCosine(k) & cosZgrid==loopCosZ(k) &...
        altitGrid==loopAltitude(k) & contamGrid==loopContam(k);
    assert(nnz(t)==size(reflMat,1),'number of matches must equal size(reflMat,1)')
    refl(t) = reflMat(:,k);
end

%return result
C.atmosphere = atmos;
C.radius = unique(radius(:))';
C.contam = unique(contamGrid(:))';
C.altitude = unique(altitGrid(:))';
C.cosZ = unique(cosZgrid(:))';
C.cosZp = unique(cosZpgrid(:))';
C.reflectance = refl;
if nargout>1
    F = griddedInterpolant({C.radius,C.contam,C.altitude,C.cosZ,C.cosZp},...
        C.reflectance,'makima','nearest');
    varargout{1} = F;
end

end