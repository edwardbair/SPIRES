function [R,varargout ] = SnowCloudSpectralRefl(varargin)
%spectral reflectance and transmittance of snow
% usage 1: R=SnowCloudSpectralRefl(substance,prop1,val1,pro2,val2,...)
% usage 2: R=SnowCloudSpectralRefl(substance,prescription,prop1,val1,pro2,val2,...)
% usage 3: R=SnowCloudSpectralRefl(prescription,prop1,val1,pro2,val2,...)
% in all cases, you can call as [R,P] = to retrieve the prescription used,
% as modified
%
%Inputs
%Substance and prescription argument, if specified, must come before the
%name-value pairs if either or both are specified
%   substance, either 'snow', 'iceCloud', 'waterCloud', or 'mixedCloud'
%       (but any unambiguous abbreviation beginning with first letter works)
%   prescription, if you want to modify an existing prescription
%
%Other inputs, name-value pairs, SEE SetSnowCloud FOR ALL POSSIBILITIES
%for example
% 'cosZ' - cosine of illumination angle, scalar, vector, or matrix
% 'radius' - effective optical radius of scatterer, scalar, vector, or matrix
%   (if not scalars, 'cosZ' and 'radius' must be same size as wavelength)
% 'wavelength' - wavelengths for which to calculate values
%
%Output
% R - structure with snow or cloud reflectance and transmittance, dimensionless,
%   same size dimensions as number of wavelengths, or radius x number of wavelengths
% P - optional, input prescription

narginchk(1,Inf)
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
iStruct = SetSnowCloud(mfilename,varargin{:});
switch iStruct.substance
    case {mixedCloud,iceCloud}
        % Mie calculations
        if iStruct.lookup % Mie parameters by lookup table
            M = lookupMie('ice',iStruct.iceRadius(:),iStruct.sizeUnit,...
                iStruct.wavelength(:),iStruct.waveUnit);
        else
            M = MieSphere(iStruct.iceRadius(:),iStruct.sizeUnit,iStruct.wavelength(:),...
                iStruct.waveUnit,'substance','ice');
        end
        switch iStruct.substance
            case mixedCloud
                % Mie calculations
                if iStruct.lookup % Mie parameters by lookup table
                    D = lookupMie('water',iStruct.waterRadius(:),iStruct.sizeUnit,...
                        iStruct.wavelength(:),iStruct.waveUnit);
                else
                    D(1,1) = MieSphere(iStruct.waterRadius(:),iStruct.sizeUnit,iStruct.wavelength(:),...
                        iStruct.waveUnit,'substance','water');
                end
                density = [iceDensity waterDensity];
                [ri,rw] = checkSizes(iStruct.iceRadius(:),iStruct.waterRadius(:));
                AllRad = [ri rw];
                contamConc(1,1) = iStruct.wetness;
            case iceCloud
                density = iceDensity;
                contamConc = [];
                AllRad = iStruct.iceRadius(:);
        end
    case snow
        if iStruct.lookup
            if iStruct.wetSnow
                M = lookupMie('wetsnow',iStruct.iceRadius(:),iStruct.sizeUnit,...
                    iStruct.wavelength(:),iStruct.waveUnit,iStruct.wetness);
            else
                M = lookupMie('snow',iStruct.iceRadius(:),iStruct.sizeUnit,...
                    iStruct.wavelength(:),iStruct.waveUnit);
            end
        else
            if iStruct.wetSnow
                CRefIn = RefractiveIndex(iStruct.wavelength(:),'ice',iStruct.waveUnit).*(1-iStruct.wetness) +...
                    RefractiveIndex(iStruct.wavelength(:),'water',iStruct.waveUnit).*iStruct.waterConc;
                M = MieSphere(iStruct.iceRadius(:),iStruct.sizeUnit,iStruct.wavelength(:),...
                    iStruct.waveUnit,'refindex',CRefIn);
            else
                M = MieSphere(iStruct.iceRadius(:),iStruct.sizeUnit,iStruct.wavelength(:),...
                    iStruct.waveUnit,'substance','ice');
            end
        end
        if iStruct.wetSnow
            density = iceDensity*(1-iStruct.wetness)+waterDensity*iStruct.wetness;
        else
            density = iceDensity;
        end
        AllRad = iStruct.iceRadius(:);
        contamConc = [];
    case waterCloud
        if iStruct.lookup
            M = lookupMie('water',iStruct.waterRadius(:),iStruct.sizeUnit,...
                iStruct.wavelength(:),iStruct.waveUnit);
        else
            M = MieSphere(iStruct.waterRadius(:),iStruct.sizeUnit,iStruct.wavelength(:),...
                iStruct.waveUnit,'substance','water');
        end
        density = waterDensity;
        contamConc = [];
        AllRad = iStruct.waterRadius(:);
end
% incorporate dust and/or soot
if iStruct.dustySnowCloud
    if iStruct.lookup
        thisD = lookupMie('dust',iStruct.dustRadius,iStruct.sizeUnit,...
            iStruct.wavelength(:),iStruct.waveUnit);
    else
        thisD = MieSphere(iStruct.dustRadius,iStruct.sizeUnit,iStruct.wavelength(:),...
            iStruct.waveUnit,'substance','dust');
    end
    if exist('D','var')
        D = cat(2,D,thisD);
    else
        D = thisD;
    end
    density = cat(2,density,dustDensity);
    contamConc = cat(2,contamConc,iStruct.dust);
    [tempR,ar] = checkSizes(iStruct.dustRadius,AllRad);
    if size(tempR,2)>1
        tempX = unique(tempR','rows');
        tempR = tempX';
    end
    AllRad = cat(2,ar,tempR);
end
if iStruct.sootySnowCloud
    if iStruct.lookup
        thisD = lookupMie('soot',iStruct.sootRadius,iStruct.sizeUnit,...
            iStruct.wavelength(:),iStruct.waveUnit);
    else
        thisD = MieSphere(iStruct.sootRadius,iStruct.sizeUnit,iStruct.wavelength(:),...
            iStruct.waveUnit,'substance','soot');
    end
    if exist('D','var')
        D = cat(2,D,thisD);
    else
        D = thisD;
    end
    density = cat(2,density,sootDensity);
    contamConc = cat(2,contamConc,iStruct.soot);
    [tempR,ar] = checkSizes(iStruct.sootRadius,AllRad);
    if size(tempR,2)>1
        tempX = unique(tempR','rows');
        tempR = tempX';
    end
    AllRad = cat(2,ar,tempR);
end

%% averaged Mie parameters if dirty snow, or mixed phase cloud
if iStruct.substance==mixedCloud || iStruct.dustySnowCloud || iStruct.sootySnowCloud
    allM = cell(1,length(M)+length(D));
    allM{1} = M;
    n = 2;
    for k=1:length(D)
        allM{n} = D(k);
        n = n+1;
    end
    conc = [1-sum(contamConc,2) contamConc];
    mix = MieMixture(allM,AllRad,density,conc);
    M = mix;
end

%% interpolate some because of incorrect Mie values
badMie = M.omega>1 | M.omega<0 | M.g>1 | M.g<0;
if nnz(badMie)
    fn = fieldnames(M);
    for k=1:length(fn)
        F = fit(iStruct.wavelength(~badMie),M.(fn{k})(~badMie),'pchipinterp');
        M.(fn{k})(badMie) = F(iStruct.wavelength(badMie));
    end
end

%% vector/matrix sizes of mu0 and omega must be equal
if isempty(iStruct.cosZ)
    mu0 = [];
    omega = M.omega;
else
    [mu0,omega] = checkSizes(iStruct.cosZ(:),M.omega);
end

if all(isinf(iStruct.WE))
    tau = Inf;
else
    switch iStruct.substance
        case mixedCloud
            mfrac = [1-iStruct.wetness iStruct.wetness];
            vfrac = mfrac./[iceDensity waterDensity];
            vfrac = vfrac/sum(vfrac);
            [ri,rw,WE,Qext,wQext] =...
                checkSizes(convertLengthUnits(iStruct.iceRadius(:),iStruct.sizeUnit,'m'),...
                convertLengthUnits(iStruct.waterRadius(:),iStruct.sizeUnit,'m'),...
                iStruct.WE,M.Qext,D(1).Qext);
            tau = tauSnow(ri,WE,Qext)*vfrac(1)+tauCloud(rw,WE,wQext)*vfrac(2);
        case {snow,iceCloud}
            [r,WE,Qext] = checkSizes(convertLengthUnits(iStruct.iceRadius(:),iStruct.sizeUnit,'m'),...
                iStruct.WE,M.Qext);
            tau = tauSnow(r,WE,Qext);
        case waterCloud
            [r,WE,Qext] =...
                checkSizes(convertLengthUnits(iStruct.waterRadius(:),iStruct.sizeUnit,'m'),...
                iStruct.WE,M.Qext);
            tau = tauCloud(r,WE,Qext);
        otherwise
            error('''substance'' ''%s'' not recognized',char(iStruct.substance))
    end
end

if isscalar(iStruct.R0)
    [g,R0,~] = checkSizes(M.g,iStruct.R0,omega);
elseif size(iStruct.R0,1)==1
    [g,R0,~] = checkSizes(M.g,mean(iStruct.R0,2),omega);
elseif istable(iStruct.R0)
    FR0 = griddedInterpolant(iStruct.R0.wavelength,mean(iStruct.R0.reflectance,2),...
        'pchip','nearest');
    [g,R0,~] = checkSizes(M.g,FR0(iStruct.wavelength(:)),omega);
else
    FR0 = griddedInterpolant(iStruct.wavelength(:),mean(iStruct.R0,2),...
        'pchip','nearest');
    [g,R0,~] = checkSizes(M.g,FR0(iStruct.wavelength(:)),omega);
end
if contains(iStruct.method,'disort')
    error('''method'' ''disort'' not yet implemented')
else
    [refl,trans,beam] = twostream(mu0,omega,g,iStruct.method,'tau',tau,'R0',R0);
end


%% fractional coverage, if applicable
if iStruct.substance==snow && ~isempty(iStruct.fractionalCoverage)
    assert(length(iStruct.fractionalCoverage)>=2,...
        'fractionalCoverage must have 2 or more values, corresponding to fSCA and 1 or more endmembers')
    if istable(iStruct.R0)
        assert(length(iStruct.fractionalCoverage) == 1+size(iStruct.R0.reflectance,2),...
            'if endmembers besides snow, need spectrum for each')
    else
        assert(length(iStruct.fractionalCoverage) == 1+size(iStruct.R0,2),...
            'if 1 endmember besides snow, neeed 1 spectrum or constant, if 2, need 2')
    end
    % make sure fractions sum to 1.0
    s = sum(iStruct.fractionalCoverage);
    if s~=1
        iStruct.fractionalCoverage = iStruct.fractionalCoverage/s;
    end
    % R0 same size as wavelength, or scalar, or interpolate, x # of non-snow endmembers
    if istable(iStruct.R0)
        for k=1:size(iStruct.R0.reflectance,2)
            FR0 = griddedInterpolant(iStruct.R0.wavelength,...
                iStruct.R0.reflectance(:,k),'pchip','nearest');
            [~,f] = checkSizes(refl,FR0(iStruct.wavelength));
            rf(:,k) = f; %#ok<AGROW>
        end
    else
        for k=1:size(iStruct.R0,2)
            [~,f] = checkSizes(refl,iStruct.R0(:,k));
            rf(:,k) = f; %#ok<AGROW>
        end
    end
    mixR = refl*iStruct.fractionalCoverage(1);
    for k=2:length(iStruct.fractionalCoverage)
        mixR = mixR+rf(:,k-1)*iStruct.fractionalCoverage(k);
    end
    refl = mixR;
end

%% results into output structure
R.refl = refl;
if any(trans(:)>0)
    R.trans = trans;
    if any(beam(:)>0)
        R.beam = beam;
    end
end
if nargout>1
    fn = fieldnames(iStruct);
    for k=1:length(fn)
        if isnumeric(iStruct.(fn{k}))
            if isscalar(unique(iStruct.(fn{k})(:)))
                iStruct.(fn{k}) = unique(iStruct.(fn{k})(:));
            end
        end
    end
    varargout{1} = iStruct;
end

end