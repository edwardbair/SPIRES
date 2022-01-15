function [I,files ] = prepBandAlbedoLookupClean(folder,substance,varargin)
% [I,files ] = prepBandAlbedoLookupClean(folder,substance,Name/Value pairs)
%create structures of snow or cloud reflectances across
%wavelength intervals or sensor bandpasses
%Generally this function is used to as input to create interpolations
%that provide rapid calculations of band-integrated reflectances
%
%Note 5/14/2020, split into 2 parts for clean and dirty snow. The dirty
%snow calculation is revised to estimate the difference between clean and
%dirty.
%
%Input
% folder - location to write output
% substance - 'snow','iceCloud','waterCloud','mixedCloud'
% Optional input, name/value pairs
% 'RefAtmos' - reference atmosphere for spectrum, default 'MLW'
% Wavelength properties, either 'bandPass' or 'sensor'/'band' must be
%   specified, but not both
% 'waveUnit' – units for bandPass or sensor/band, default 'nm'
% 'bandPass' – Nx2 matrix for multispectral sensors, or use [.28 4] or
%   [280 4000] to get full spectrum albedo
% 'sensor' – instead of 'wavelength', can specify any multispectral sensor
%   in the SensorTable.m function
% 'bands' – bands of the sensor, either numeric vector, cell vector, or
%   categorical vector of bands, or, if omitted, all bands for that sensor
% 'nRadii' - number of radii over which to calculate over the limits of
%   SnowCloudLimits, default 35 for clean snow
% 'mu0' and 'muS' - cosines of illumination angle on flat and on slopes,
%   default [.05:.025:.25 .3:.05:1] for flat, zero included for slopes
% 'elevation' - vector of altitudes (in meters) for which to calculate
%   values, default [2000 3000 4000 5500 6000 7000], weirdly spaced because
%   of change in atmospheric attenuation around 5500 m
% (will vary because the spectral distribution of irradiance varies with
% altitude)
%
%Output
% I - structure containing the inputs
% file - where results are written
%
%Note: these calculations are for open surfaces on slopes. The spectrum of
%diffuse irradiance in the trees would have more in the near-infrared
%because of scattering from the leaves or needles.

%% parse input and check sizes
narginchk(1,17)
nargoutchk(0,2)
bandValidation = @(x) ((isnumeric(x) && all(x(:)>0)) || iscell(x) || iscategorical(x));
bpValidation = @(x) isnumeric(x) && all(x(:)>=0);
p = inputParser;
addRequired(p,'folder',@ischar)
addRequired(p,'substance',@ischar)
addParameter(p,validatestring('bandpass',{'bandp','bandpass'}),...
    [],bpValidation)
addParameter(p,validatestring('waveunit',{'waveu','waveunit','lambda'}),...
    'nm',@ischar)
addParameter(p,validatestring('sensor',{'sens','sensor'}),...
    '',@(x) ischar(x) || iscategorical(x))
addParameter(p,validatestring('bands',{'band','bands'}),...
    [],bandValidation)
addParameter(p,'refatmos','mlw',@ischar)
addParameter(p,'elevation',...
    [2000 3000 4000 5500 6000 7000],@isnumeric)
addParameter(p,'cosZ',[.05:.025:.25 .3:.05:1],@isnumeric)
addParameter(p,'nradii',35,@isnumeric)
parse(p,folder,substance,varargin{:});

SnowStruct = SnowCloudLimits;
folder = p.Results.folder;

%% just one atmosphere available now
assert(strcmpi(p.Results.refatmos,'mlw'),'just ''mlw'' atmosphere available now')
waveUnit = p.Results.waveunit;

%% bandpass or sensor?
assert(xor(isempty(p.Results.bandpass),isempty(p.Results.sensor)),...
    'must specify ''bandPass'' or ''sensor'', but not both')
if isempty(p.Results.bandpass)
    sensT = SensorTable(p.Results.sensor,waveUnit);
    if isempty(p.Results.bands)
        bandPass = [sensT.LowerWavelength sensT.UpperWavelength];
        bands = sensT.Band;
    else
        bTbl = table(categorical(p.Results.bands(:)),'VariableNames',{'Band'});
        nTbl = innerjoin(sensT,bTbl);
        bandPass = [nTbl.LowerWavelength nTbl.UpperWavelength];
        bands = nTbl.Band;
    end
else
    bandPass = p.Results.bandpass;
end

% lookup table uses nm, so convert
bandPass = convertLengthUnits(bandPass,waveUnit,'nm');

% values of mu0, muS
ucz = p.Results.cosZ;
ucS = [0 ucz];

%% different radius limits depending on substance
if contains(p.Results.substance,'snow','IgnoreCase',true)
    substance = 'snow';
    radii = SnowStruct.snowRadius;
elseif contains(p.Results.substance,'icec','IgnoreCase',true)
    radii = SnowStruct.iceCloudRadius;
    substance = 'iceCloud';
elseif contains(p.Results.substance,'wate','IgnoreCase',true)
    radii = SnowStruct.waterCloudRadius;
    substance = 'waterCloud';
elseif contains(p.Results.substance,'mixe','IgnoreCase',true)
    substance = 'mixedCloud';
    radii = SnowStruct.iceCloudRadius;
else
    error('substance ''%s'' not recognized',p.Results.substance)
end
ucr = linspace(sqrt(radii(1)),sqrt(radii(2)),p.Results.nradii).^2;
ucssa = sort(radius2SSA(ucr,'um'));

%% inputs into the I structure
I.elevation = p.Results.elevation;
I.mu0 = ucz;
I.SSA = ucssa;
I.radius = ucr;
I.muS = ucS;
I.bandPass = bandPass;
if exist('bands','var')
    I.bands = bands;
end

%% reflectance wavelengths use those at which refractive index is measured
if strcmpi(substance,'water')
    [~,wavelength] = RefractiveIndex([],'water','nm');
else
    [~,wavelength] = RefractiveIndex([],'ice','nm');
end
wavelength = wavelengthsNeeded(wavelength,bandPass);

%solar wavelengths available for this set of bandPasses based on SMARTS
altit = p.Results.elevation;
SP = defaultSMARTSinput(p.Results.refatmos,'altit',altit(1),'cosZ',ucz(end));
S = SMARTSMain(SP);
% solar wavelengths needed
waveL = wavelengthsNeeded(S.spectralTbl.waveL,...
    [min(wavelength) max(wavelength)]);
tmpT = table(waveL);
STbl = innerjoin(S.spectralTbl,tmpT);
wlmn = STbl.waveL(1);
wlmx = STbl.waveL(end);

%N-D grids of the inputs for clean snow
[V.elevation,V.mu0,V.SSA,V.muS] = ndgrid(I.elevation,I.mu0,I.SSA,I.muS);

%% main nested loop
fprintf('size of loop %d -- %d elevations, %d solar angles\n',...
    numel(V.elevation),length(I.elevation),length(I.mu0));
for ka=1:length(I.elevation)
    for kz=1:length(I.mu0)
        % irradiance depends on elevation and cosZ
        SP = defaultSMARTSinput(p.Results.refatmos,'altit',I.elevation(ka),...
            'cosZ',I.mu0(kz),'WLMN',wlmn,'WLMX',wlmx);
        smartsS = SMARTSMain(SP);
        
        %table of solar direct and diffuse irradiance at this elevation and
        %solar zenith angle, covering just the wavelengths of the bandPasses
        solarTbl = table(smartsS.spectralTbl.waveL,...
            [smartsS.spectralTbl.HorzDirect smartsS.spectralTbl.HorzDiffuse],...
            'VariableNames',{'wavelength','irradiance'});
        solarTbl.Properties.VariableUnits = {'nm','W m^(-2) nm^(-1)'};
        
        %trigger a call to bandPassReflectance to set the solar interpolation
        %functions (okay to use random numbers for the reflectances)
        reflTbl = table(wavelength,rand(length(wavelength),2),...
            'VariableNames',{'wavelength','reflectance'});
        reflTbl.Properties.VariableUnits = {'nm',''};
        [~,pscript] = bandPassReflectance(reflTbl,solarTbl,'bandPass',...
            bandPass,'mu0',I.mu0(1),'muS',I.muS(1));
        
        %% loop through combination
        t = V.elevation==I.elevation(ka) & V.mu0==I.mu0(kz);
        subNames = {'SSA','muS'};
        for v=1:length(subNames)
            X = V.(subNames{v})(t);
            subV.(subNames{v}) = X(:);
        end
        % pass subV to the workers
        passVals = {subV};
        W = parallel.pool.Constant(passVals);
        nLoop = numel(subV.SSA);
        % allocate temporary reflectance array
        %transpose output array because MATLAB use column-major order
        reflArray = zeros(numel(subV.SSA),size(bandPass,1)).';
        loop_mu0 = I.mu0(kz);
        % this loop computes the clean-snow values
        parfor kr=1:nLoop
            pv = W.Value;
            % diffuse reflectance is independent of illumination angle
            P = fwdPrescription(substance,'cosZ',[],'SSA',...
                pv{1}.SSA(kr),'lookup',true,'sizeU','um',...
                'waveU','nm','wavelength',wavelength);
            [~,Rdif] =  SPIReS_fwd(P);
            %direct reflectance depends on local illumination angle
            P = fwdPrescription(substance,'cosZ',pv{1}.muS(kr),'SSA',...
                pv{1}.SSA(kr),'lookup',true,'sizeU','um',...
                'waveU','nm','wavelength',wavelength);
            [~,Rb] =  SPIReS_fwd(P);
            reflTbl = table(wavelength,[Rb.reflectance Rdif.reflectance],...
                'VariableNames',{'wavelength','reflectance'});
            reflTbl.Properties.VariableUnits = {'nm',''};
            bandRefl = bandPassReflectance(reflTbl,pscript,...
                'mu0',loop_mu0,'muS',pv{1}.muS(kr));
            reflArray(:,kr) = bandRefl.'; %#ok<PFOUS>
        end
        % save as a checkpoint
        checkname = sprintf('checkClean_%d_%d.mat',ka,kz);
        save(fullfile(folder,checkname),'ka','kz','reflArray','-v7.3')
        if mod(ka*kz,10)==1
            fprintf('loop elevation %d angle %d completed\n',ka,kz);
        end
    end
end
files = dir(fullfile(folder,'checkClean*.mat'));
end