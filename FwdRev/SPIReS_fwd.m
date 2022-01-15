function [refl,varargout] = SPIReS_fwd(prescription)
%spectral reflectance and transmittance of snow
% refl = SPIReS_fwd(prescription)
% [refl,Tbl] = SPIReS_fwd(prescription)
%   (use setPrescription to set prescription)
% If the prescription.RadiativeTransfer.diffuseFraction function is set,
% modeled reflectance is corrected for terrain and diffuse illumination
%
%Output
% refl - spectral reflectance as a vector corresponding to wavelengths
%Optional output
% Tbl - table with snow or cloud reflectance and transmittance, dimensionless,
%   same size dimensions as number of wavelengths, or radius x number of wavelengths
%   as set in the prescription

narginchk(1,1)
nargoutchk(0,2)

%% all substances
% all categories
snow = categorical({'snow'});
iceCloud = categorical({'iceCloud'});
waterCloud = categorical({'waterCloud'});
mixedCloud = categorical({'mixedCloud'});

%% parse inputs
iceDensity = 917; % kg/m^3
waterDensity = 1000;
dustDensity = 2800;
sootDensity = 1800;
Illumination = prescription.Illumination;
Spectrum = prescription.Spectrum;
Substrate = prescription.Substrate;
RadiativeTransfer = prescription.RadiativeTransfer;
switch prescription.substance
    case {mixedCloud,iceCloud}
        cloudStruct = char(prescription.substance);
        % Mie calculations
        if RadiativeTransfer.lookupMie % Mie parameters by lookup table
            M = lookupMie('ice',prescription.(cloudStruct).radius,...
                prescription.(cloudStruct).sizeUnit,...
                Spectrum.wavelength,Spectrum.waveUnit);
        else
            M = mieSphere(prescription.(cloudStruct).radius,...
                prescription.(cloudStruct).sizeUnit,Spectrum.wavelength,...
                Spectrum.waveUnit,'substance',RadiativeTransfer.icer);
        end
        switch prescription.substance
            case mixedCloud
                % Mie calculations
                if RadiativeTransfer.lookupMie % Mie parameters by lookup table
                    D = lookupMie('water',prescription.mixedCloud.waterRadius,...
                        prescription.mixedCloud.sizeUnit,...
                        Spectrum.wavelength,Spectrum.waveUnit);
                else
                    D(1,1) = mieSphere(prescription.mixedCloud.waterRadius,...
                        prescription.mixedCloud.sizeUnit,Spectrum.wavelength,...
                        Spectrum.waveUnit,'substance','water');
                end
                density = [iceDensity waterDensity];
                [ri,rw] = checkSizes(prescription.mixedCloud.radius,prescription.mixedCloud.waterRadius);
                AllRadii = [ri rw];
                contamConc(1,1) = prescription.mixedCloud.wetness;
            case iceCloud
                density = iceDensity;
                contamConc = [];
                AllRadii = prescription.iceCloud.radius;
        end
    case snow
        if RadiativeTransfer.lookupMie
            if prescription.snow.wetSnow
                M = lookupMie('wetsnow',prescription.snow.radius,prescription.snow.sizeUnit,...
                    Spectrum.wavelength,prescription.Spectrum.waveUnit,prescription.snow.wetness);
            else
                M = lookupMie('snow',prescription.snow.radius,prescription.snow.sizeUnit,...
                    Spectrum.wavelength,Spectrum.waveUnit);
            end
        else
            if prescription.snow.wetSnow
                CRefIn = RefractiveIndex(Spectrum.wavelength,...
                    RadiativeTransfer.icer,Spectrum.waveUnit).*...
                    (1-prescription.snow.wetness) +...
                    RefractiveIndex(Spectrum.wavelength,'water',...
                    Spectrum.waveUnit).*prescription.snow.wetness;
                M = mieSphere(prescription.snow.radius,prescription.snow.sizeUnit,...
                    Spectrum.wavelength,...
                    Spectrum.waveUnit,'refindex',CRefIn);
            else
                M = mieSphere(prescription.snow.radius,prescription.snow.sizeUnit,...
                    Spectrum.wavelength,Spectrum.waveUnit,...
                    'substance',RadiativeTransfer.icer);
            end
        end
        if prescription.snow.wetSnow
            density = iceDensity*(1-prescription.snow.wetness)+...
                waterDensity*prescription.snow.wetness;
        else
            density = iceDensity;
        end
        AllRadii = prescription.snow.radius;
        contamConc = [];
    case waterCloud
        if RadiativeTransfer.lookupMie
            M = lookupMie('water',prescription.waterCloud.radius,...
                prescription.waterCloud.sizeUnit,...
                Spectrum.wavelength,Spectrum.waveUnit);
        else
            M = mieSphere(prescription.watewaterCloud.radius,...
                prescription.waterCloud.sizeUnit,Spectrum.wavelength,...
                Spectrum.waveUnit,'substance','water');
        end
        density = waterDensity;
        contamConc = [];
        AllRadii = prescription.waterCloud.radius;
end
% incorporate dust and/or soot
if prescription.substance==snow && ~prescription.snow.cleanSnow
    for k=1:prescription.snow.nLAP
        if RadiativeTransfer.lookupMie
            thisD = lookupMie(prescription.snow.LAP{k},prescription.snow.LAPradius(k),...
                prescription.snow.sizeUnit,...
                Spectrum.wavelength,Spectrum.waveUnit);
        elseif iscell(prescription.snow.LAP)
            thisD = mieSphere(prescription.snow.LAPradius(k),...
                prescription.snow.sizeUnit,Spectrum.wavelength,...
                Spectrum.waveUnit,'substance',prescription.snow.LAP{k});
        elseif isnumeric(prescription.snow.LAP)
            thisD = mieSphere(prescription.snow.LAPradius(k),...
                prescription.snow.sizeUnit,Spectrum.wavelength,...
                Spectrum.waveUnit,'refindex',prescription.snow.LAP(:,k));
        else
            error('class of prescription.snow.LAP (%s) not recognized',...
                class(prescription.snow.LAP))
        end
        if exist('D','var')
            D = cat(2,D,thisD);
        else
            D = thisD;
        end
        if isnumeric(prescription.snow.LAP) || strcmpi(prescription.snow.LAP{k},'dust')
            density = cat(2,density,dustDensity);
        else
            density = cat(2,density,sootDensity);
        end
        contamConc = cat(2,contamConc,prescription.snow.LAPfraction(k));
        [tempR,ar] = checkSizes(prescription.snow.LAPradius(k),AllRadii);
        if size(tempR,2)>1
            tempX = unique(tempR','rows');
            tempR = tempX';
        end
        AllRadii = cat(2,ar,tempR);
    end
end

%% averaged Mie parameters if dirty snow, or mixed phase cloud
if prescription.substance==mixedCloud || (prescription.substance==snow && ~prescription.snow.cleanSnow)
    allM = cell(1,length(M)+length(D));
    allM{1} = M;
    n = 2;
    for k=1:length(D)
        allM{n} = D(k);
        n = n+1;
    end
    conc = [1-sum(contamConc,2) contamConc];
    mix = mieMixture(allM,AllRadii,density,conc);
    M = mix;
end

%% vector/matrix sizes of mu0 and omega must be equal
if isempty(prescription.Illumination.muS)
    mu0 = [];
    omega = M.omega;
else
    [mu0,omega] = checkSizes(prescription.Illumination.muS,M.omega);
end

if prescription.substance==snow && prescription.snow.deepSnow
    tau = Inf;
else
    switch prescription.substance
        case mixedCloud
            mfrac = [1-prescription.mixedCloud.wetness prescription.mixedCloud.wetness];
            vfrac = mfrac./[iceDensity waterDensity];
            vfrac = vfrac/sum(vfrac);
            [ri,rw,WE,Qext,wQext] =...
                checkSizes(convertLengthUnits(prescription.mixedCloud.radius,prescription.mixedCloud.sizeUnit,'m'),...
                convertLengthUnits(prescription.mixedCloud.waterRadius,prescription.mixedCloud.sizeUnit,'m'),...
                prescription.mixedCloud.waterEquivalent,M.Qext,D(1).Qext);
            tau = tauSnow(ri,WE,Qext)*vfrac(1)+tauCloud(rw,WE,wQext)*vfrac(2);
        case {snow,iceCloud}
            substruct = char(prescription.substance);
            [r,WE,Qext] = checkSizes(convertLengthUnits(prescription.(substruct).radius,...
                prescription.(substruct).sizeUnit,'m'),prescription.(substruct).waterEquivalent,M.Qext);
            tau = tauSnow(r,WE,Qext);
        case waterCloud
            [r,WE,Qext] =...
                checkSizes(convertLengthUnits(prescription.waterCloud.radius,prescription.waterCloud.sizeUnit,'m'),...
                prescription.waterCloud.waterEquivalent,M.Qext);
            tau = tauCloud(r,WE,Qext);
        otherwise
            error('''substance'' ''%s'' not recognized',char(prescription.substance))
    end
end

if contains(RadiativeTransfer.method,'disort')
    error('''method'' ''disort'' not yet implemented')
end
if ~isinf(tau)
    if Substrate.useFunction
        thisR0 = zeros(length(Spectrum.wavelength),length(Substrate.fR0));
        for k=1:size(thisR0,2)
            F = Substrate.fR0{k};
            thisR0(:,k) = F(Spectrum.wavelength);
        end
        [g,R0,~] = checkSizes(M.g,mean(thisR0,2),omega);
    else
        [g,R0,~] = Substrate.R0;
    end
    [refl,trans,beam] = twostream(mu0,omega,g,RadiativeTransfer.method,'tau',tau,'R0',R0);
else
    [g,~] = checkSizes(M.g,omega);
    [refl,trans,beam] = twostream(mu0,omega,g,RadiativeTransfer.method);
end

% if fSCA<1, linear mixture
if prescription.substance==snow &&...
        (isfield(prescription.snow,'fSCA') && prescription.snow.fSCA(1)~=1)
    thisR0 = zeros(length(Spectrum.wavelength),length(Substrate.fR0));
    for k=1:size(thisR0,2)
        F = Substrate.fR0{k};
        thisR0(:,k) = F(Spectrum.wavelength);
    end
    refl = linearMixture([refl thisR0],prescription.snow.fSCA);
end

% if diffuse fraction indicated, correct for terrain
if prescription.RadiativeTransfer.calcTerrain
    refl = terrainCorrection(prescription,refl);
end

%% results into output table
if nargout>1
    if isempty(Illumination.cosZ)
        if any(contains(fieldnames(Spectrum),'waveBand'))
            R = table(Spectrum.wavelength,Spectrum.waveBand,refl,double.empty([length(Spectrum.wavelength) 0]),...
                'VariableNames',{'wavelength','waveBand''reflectance','cosZ'});
            R.Properties.VariableUnits = {Spectrum.waveUnit,'','',''};
        else
            R = table(Spectrum.wavelength,refl,double.empty([length(Spectrum.wavelength) 0]),...
                'VariableNames',{'wavelength','reflectance','cosZ'});
            R.Properties.VariableUnits = {Spectrum.waveUnit,'',''};
        end
    elseif Illumination.cosZ==Illumination.muS
        if any(contains(fieldnames(Spectrum),'waveBand'))
            [cz,w,b] = checkSizes(Illumination.cosZ,Spectrum.wavelength,Spectrum.waveBand);
            R = table(w,b,refl,cz,...
                'VariableNames',{'wavelength','waveBand','reflectance','cosZ'});
            R.Properties.VariableUnits = {Spectrum.waveUnit,'','',''};
        else
            [cz,w] = checkSizes(Illumination.cosZ,Spectrum.wavelength);
            R = table(w,refl,cz,...
                'VariableNames',{'wavelength','reflectance','cosZ'});
            R.Properties.VariableUnits = {Spectrum.waveUnit,'',''};
        end
    else
        if any(contains(fieldnames(Spectrum),'waveBand'))
            [cz,w,b] = checkSizes(Illumination.cosZ,Spectrum.wavelength,Spectrum.waveBand);
            muS = ones(size(cz))*Illumination.muS;
            R = table(w,b,refl,cz,muS,...
                'VariableNames',{'wavelength','waveBand','reflectance','cosZ','muS'});
            R.Properties.VariableUnits = {Spectrum.waveUnit,'','','',''};
        else
            [cz,w] = checkSizes(Illumination.cosZ,Spectrum.wavelength);
            muS = ones(size(cz))*Illumination.muS;
            R = table(w,refl,cz,muS,...
                'VariableNames',{'wavelength','reflectance','cosZ','muS'});
            R.Properties.VariableUnits = {Spectrum.waveUnit,'','',''};
        end
    end
    subStruct = char(prescription.substance);
    [rr,~] = checkSizes(prescription.(subStruct).radius,R.wavelength);
    newTbl = table(rr,'VariableNames',{'radius'});
    newTbl.Properties.VariableUnits = {prescription.(subStruct).sizeUnit};
    R = [R newTbl];
    if any(trans>0)
        newTbl = table(trans,'VariableNames',{'transmittance'});
        newTbl.Properties.VariableUnits = {''};
        R = [R newTbl];
        if any(beam>0)
            newTbl = table(beam,'VariableNames',{'beamTrans'});
            newTbl.Properties.VariableUnits = {''};
            R = [R newTbl];
        end
        newTbl = table(repmat(prescription.(subStruct).waterEquivalent,height(R),1),...
            'VariableNames',{'waterEquivalent'});
        newTbl.Properties.VariableUnits = {'mm'};
        R = [R newTbl];
        newTbl = table(thisR0,'VariableNames',{'R0'});
        newTbl.Properties.VariableUnits = {''};
        R = [R newTbl];
    elseif prescription.substance==snow && prescription.snow.fSCA(1)<1
        newTbl = table(thisR0,repmat(prescription.snow.fSCA,height(R),1),...
            'VariableNames',{'R0','fSCA_and'});
        newTbl.Properties.VariableUnits = {'',''};
        R = [R newTbl];
    end
    if any(contains(R.Properties.VariableNames,'waveBand'))
        R = sortrows(R,{'cosZ','radius','waveBand','wavelength'});
    else
        R = sortrows(R,{'cosZ','radius','wavelength'});
    end
    varargout{1} = R;
end
end