function [ T,varargout] = SnowCloudIntgRefl(varargin)
% [ T [,P]] = SnowCloudIntgRefl(solarTbl,substance,[prescription and/or prop/value pairs\)
% [ T [,P]] = SnowCloudIntgRefl(substance,[prescription and/or prop/value pairs\)
%band-integrated reflectance and transmittance of snow or cloud, coupling
%spectral models of snow/cloud reflectance and, optionally, incoming solarradiation
%in the calculation of snow reflectance across specified wavelenth bands
%
%Optional input
% solarTbl
%   If the wavelength range(s) is/are broad enough that there is spectral
%   variability of the irradiance across the band(s), then solar
%   irradiance must be accounted for, in a table
%   The first column in the table is wavelength, for the irradiance
%   spectrum
%   The second column is irradiance, either a vector for the value at each
%   wavelength, or a matrix with 2 columns where the first is the direct
%   irradiance and the next is the diffuse irradiance.
%   The Table.Properties.VariableUnits must be specified, with the first
%   some length unit for the wavelengths, typically 'mum' or 'nm' and the
%   second unit is ''
%
% If solarTbl is not the first argument, then band-averaged reflectivity is
% used.
%
%The following variable inputs specify the snow or cloud properties, by calling
%SetSnowCloud (see that function for additional details)
%
%The first values must be the substance, either 'snow', 'iceCloud',
%'waterCloud', or 'mixedCloud' (any unambiguous abbreviation beginning
%with first letter works) or you can specify the prescription output from
%SetSnowCloud. The other inputs are prop/value pairs
% Wavelength properties, either 'bandPass' or 'sensor'/'band' must be
%   specified, but not both
% 'waveUnit' – units for bandPass or sensor/band, default same as the units
%   in the solar table, otherwise no default (to prevent errors)
% 'bandPass' – Nx2 matrix for multispectral sensors, or use [.28 4] or
%   [280 4000] to get full spectrum albedo
% 'sensor' – instead of 'wavelength', can specify any sensor in the
%   SensorTable.m function
% 'bands' – bands of the sensor, either numeric vector, cell vector, or
%   categorical vector of bands,or, if omitted, all bands for that sensor
% 'R0' - reflectance of underlying surface or non-snow part of pixel,
%   scalar or table with columns wavelength and reflectance and the
%   wavelength units specified in Table.Properties.VariableUnits
%
%Values for 'radius' and 'cosZ' can be specified as scalars or vectors, but if
%vectors they must be the same size. To get all combinations of radius/cosZ
%values, use meshgrid or ndgrid first and convert the results to vectors.
%If solar Table is an argument, you must use the same solar geometry in the
%software that generates the irradiance values.
%
%Output
% T - table of snow or cloud reflectance and transmittance, dimensionless,
%   same height as number of bands x number of radii (which equals number
%   of cosZ)
%Optional output
% P - prescription of snow/cloud properties

%%
narginchk(2,Inf)
nargoutchk(0,2)

fillSolarWaveUnits = false;
if istable(varargin{1})
    solarTbl = varargin{1};
    newVarargin = varargin(2:end);
    % add default wavelength if not specified
    waveSpecified = false;
    for k=1:length(newVarargin)
        if ischar(newVarargin{k})
            waveSpecified = contains(newVarargin{k},'waveu','IgnoreCase',true);
        end
    end
    if waveSpecified
        if isempty(solarTbl.Properties.VariableUnits)
            fillSolarWaveUnits = true;
        end
    else
        solarWaveUnit = solarTbl.Properties.VariableUnits{1};
        assert(~isempty(solarWaveUnit),...
            'input solar table must have wavelength units specified in the Properties.VariableUnits')
        newVarargin{end+1} = 'waveUnit';
        newVarargin{end+1} = solarWaveUnit;
    end
else
    newVarargin = varargin;
end

SnowCloudP = SetSnowCloud(mfilename,newVarargin{:});
assert(~SnowCloudP.spectrometer,...
    'code %s designed for spectrally integrated measurements, not a spectrometer',...
    mfilename);
if fillSolarWaveUnits
    solarTbl.Properties.VariableUnits = {SnowCloudP.waveUnit};
elseif exist('solarTbl','var')
    assert(strcmp(solarTbl.Properties.VariableUnits{1},SnowCloudP.waveUnit),...
        'if ''waveUnits'' and solarTbl.Properties.VariableUnits are both specified, they must be the same')
end

% generate reflectance spectrum for the snow or cloud
% hold the wavelength ranges for later integration
minW = min(SnowCloudP.bandPass(:));
maxW = max(SnowCloudP.bandPass(:));
if contains(char(SnowCloudP.substance),'water','IgnoreCase',true)
    [~,wv] = RefractiveIndex([],'water',SnowCloudP.waveUnit);
else
    [~,wv] = RefractiveIndex([],'ice',SnowCloudP.waveUnit);
end
k1 = find(wv<=minW,1,'last');
k2 = find(wv>=maxW,1,'first');
wv = wv(k1:k2);
P1 = SnowCloudP;

% n-d grid of radius cosZ and wavelength
snow = categorical({'snow'});
waterCloud = categorical({'waterCloud'});
iceCloud = categorical({'iceCloud'});
mixedCloud = categorical({'mixedCloud'});
switch P1.substance
    case {snow,iceCloud,mixedCloud}
        [r,cz,wave] = ndgrid(unique(P1.iceRadius),unique(P1.cosZ),wv);
        P1.iceRadius = r(:);
    case waterCloud
        [r,cz,wave] = ndgrid(unique(P1.waterRadius),unique(P1.cosZ),wv);
        P1.waterRadius = r(:);
    otherwise
        error('substance ''%s'' not recognized',char(P1.substance))
end
P1.cosZ = cz(:);
P1.wavelength = wave(:);

% direct reflectance
M = SnowCloudSpectralRefl(P1);
refl = M.refl;
% diffuse reflectance
if exist('solarTbl','var') && size(solarTbl.irradiance,2)==2
    P2 = P1;
    P2.cosZ = [];
    M = SnowCloudSpectralRefl(P2);
    refl = [refl M.refl];
end

% spectral reflectance for each radius, cosZ combination
cosZ = unique(P1.cosZ);
switch P1.substance
    case {snow,iceCloud,mixedCloud}
        rad = unique(P1.iceRadius);
        radName = 'iceRadius';
    case waterCloud
        rad = unique(P1.waterRadius);
        radName = 'waterRadius';
end
for kr=1:length(rad)
    for kc=1:length(cosZ)
        thisPair = P1.(radName)==rad(kr) & P1.cosZ==cosZ(kc);
        thisRefl = refl(thisPair,:);
        w = P1.wavelength(thisPair);
        
        for b=1:size(P1.bandPass,1)
            if exist('solarTbl','var')
                if any(thisRefl(:)<0)
                    pause
                end
                R = bandPassReflectance(w,P1.waveUnit,thisRefl,solarTbl.irradiance,...
                    'bandPass',P1.bandPass(b,:),'irradWavelength',solarTbl.wavelength);
            else
                R = bandPassReflectance(w,P1.waveUnit,thisRefl,ones(size(thisRefl)),...
                    'bandPass',P1.bandPass(b,:));
            end
            if ~isempty(P1.sensor)
                thisTbl = table(P1.sensor,P1.bands(b),P1.bandPass(b,:),...
                    rad(kr),cosZ(kr),P1.WE,R,...
                    'VariableNames',...
                    {'Sensor','Band','bandPass',radName,'cosZ','waterEquivalent','reflectance'});
                thisTbl.Properties.VariableUnits =...
                    {'','',P1.waveUnit,P1.sizeUnit,'',P1.weUnit,''};
            else
                thisTbl = table(P1.bandPass(b,:),rad(kr),cosZ(kc),P1.WE,R,'VariableNames',...
                    {'bandPass',radName,'cosZ','waterEquivalent','reflectance'});
                thisTbl.Properties.VariableUnits =...
                    {P1.waveUnit,P1.sizeUnit,'',P1.weUnit,''};
            end
            if ~isempty(P1.dust)
                addTbl = table(P1.dust,P1.dustRadius,'VariableNames',...
                    {'dustConc','dustRadius'});
                addTbl.Properties.VariableUnits = {'',P1.sizeUnit};
                thisTbl = [thisTbl addTbl]; %#ok<AGROW>
            end
            if ~isempty(P1.soot)
                addTbl = table(P1.soot,P1.sootRadius,'VariableNames',...
                    {'sootConc','sootRadius'});
                addTbl.Properties.VariableUnits = {'',P1.sizeUnit};
                thisTbl = [thisTbl addTbl]; %#ok<AGROW>
            end
            if ~isempty(P1.wetness) && P1.substance~=waterCloud
                addTbl = table(P1.wetness,unique(P1.waterRadius),'VariableNames',...
                    {'wetness','waterRadius'});
                addTbl.Properties.VariableUnits = {'',P1.sizeUnit};
                thisTbl = [thisTbl addTbl]; %#ok<AGROW>
            end
            if ~isempty(P1.fractionalCoverage)
                addTbl = table(P1.fractionalCoverage,'VariableNames',...
                    {'fractionalCoverage'});
                addTbl.Properties.VariableUnits = {''};
                thisTbl = [thisTbl addTbl]; %#ok<AGROW>
            end
            if kr==1 && kc==1 && b==1
                T = thisTbl;
            else
                T = [T; thisTbl]; %#ok<AGROW>
            end
        end
        radCosPair(n,:) = [rad(kr) cosZ(kc)];
         reflMatrix(n,:) = thisRefl(:)';
    end
end

%all bandpasses and wavelengths at once
if exist('solarTbl','var')
    if isempty(P1.sensor)
        R = bandPassRefl(w,P1.waveUnit,reflMatrix,solarTbl,...
            'bandPass',P1.bandPass,'useParallel',P1.useParallel);
    else
        R = bandPassRefl(w,P1.waveUnit,reflMatrix,solarTbl,...
            'sensor',P1.sensor,'band',P1.bands,'useParallel',P1.useParallel);
    end
else
    if isempty(P1.sensor)
        R = bandPassRefl(w,P1.waveUnit,reflMatrix,[],...
            'bandPass',P1.bandPass,'useParallel',P1.useParallel);
    else
        R = bandPassRefl(w,P1.waveUnit,reflMatrix,[],...
            'sensor',P1.sensor,'band',P1.bands,'useParallel',P1.useParallel);
    end
    
end

%results into table
rad = radCosPair(:,1);
cosZ = radCosPair(:,2);
T = table;
for k=1:length(rad)
    thisRad = R(k,:)';
    if ~isempty(P1.bands) %bands is usually empty so had to change
        theseBands = P1.bands(:);
    else
       theseBands=[1:size(P1.bandPass,1)]'; 
    end
    thisSize = repmat(rad(k),size(theseBands));
    thisCos = repmat(cosZ(k),size(thisSize));
    sens = repmat(P1.sensor,size(thisSize));
    we = repmat(P1.WE,size(thisSize));
    if ~isempty(P1.sensor)
        thisTbl = table(sens,theseBands,P1.bandPass,thisSize,...
            thisCos,we,thisRad,...
            'VariableNames',...
            {'Sensor','Band','bandPass',radName,'cosZ',...
            'waterEquivalent','reflectance'});
        thisTbl.Properties.VariableUnits =...
            {'','',P1.waveUnit,P1.sizeUnit,'',P1.weUnit,''};
    else
        thisTbl = table(P1.bandPass,thisSize,thisCos,we,thisRad,...
            'VariableNames',...
            {'bandPass',radName,'cosZ','waterEquivalent','reflectance'});
        thisTbl.Properties.VariableUnits =...
            {P1.waveUnit,P1.sizeUnit,'',P1.weUnit,''};
    end
    if ~isempty(P1.dust)
        addTbl = table(repmat(P1.dust,size(thisSize)),...
            repmat(P1.dustRadius,size(thisSize)),...
            'VariableNames',...
            {'dustConc','dustRadius'});
        addTbl.Properties.VariableUnits = {'',P1.sizeUnit};
        thisTbl = [thisTbl addTbl]; %#ok<AGROW>
    end
    if ~isempty(P1.soot)
        addTbl = table(repmat(P1.soot,size(thisSize)),...
            repmat(P1.sootRadius,size(thisSize)),...
            'VariableNames',...
            {'sootConc','sootRadius'});
        addTbl.Properties.VariableUnits = {'',P1.sizeUnit};
        thisTbl = [thisTbl addTbl]; %#ok<AGROW>
    end
    if ~isempty(P1.wetness) && P1.substance~=waterCloud
        addTbl = table(repmat(P1.wetness,size(thisSize)),...
            repmat(unique(P1.waterRadius(:)),size(thisSize)),...
            'VariableNames',...
            {'wetness','waterRadius'});
        addTbl.Properties.VariableUnits = {'',P1.sizeUnit};
        thisTbl = [thisTbl addTbl]; %#ok<AGROW>
    end
    if ~isempty(P1.fractionalCoverage)
        addTbl = table(repmat(P1.fractionalCoverage,size(thisSize)),...
            'VariableNames',...
            {'fractionalCoverage'});
        addTbl.Properties.VariableUnits = {''};
        thisTbl = [thisTbl addTbl]; %#ok<AGROW>
    end
    T = [T; thisTbl]; %#ok<AGROW>
end

if nargout>1
    varargout{1} = P1;
end
end