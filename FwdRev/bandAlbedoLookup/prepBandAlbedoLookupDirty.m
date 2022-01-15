function [I,files ] = prepBandAlbedoLookupDirty(folder,varargin)
% [I,files ] = prepBandAlbedoLookupDirty(folder,Name/Value pairs)
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
%   SnowCloudLimits, default is 12 for dirty snow
% 'mu0' and 'muS' - cosines of illumination angle on flat and on slopes,
%   default [.3:.1:.8] for flat and slopes
% 'nLAP' - number of concentrations for the light-absorbing particles,
%   default 12 for soot, 1.5x for dust
% 'elevation' - vector of altitudes (in meters) for which to calculate
%   values, default [3000 5500 6000] with 'nearest' extrapolation for dirty
%   snow
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
    [3000 5500 6000],@isnumeric)
addParameter(p,'cosZ',(.3:.1:.8),@isnumeric)
addParameter(p,'nradii',12,@isnumeric)
addParameter(p,'nlap',12,@isnumeric)
parse(p,folder,varargin{:});

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

%% radius limits
radii = SnowStruct.snowRadius;
radii(1) = radii(1)*2;
radii(2) = radii(2)/2;
ucr = linspace(sqrt(radii(1)),sqrt(radii(2)),p.Results.nradii).^2;
ucssa = sort(radius2SSA(ucr,'um'));

%% LAP characteristics
nLAP = p.Results.nlap;
soot = linspace(min(SnowStruct.soot),max(SnowStruct.soot),nLAP);
dust = linspace(min(SnowStruct.dust),max(SnowStruct.dust),round(1.3*nLAP));
sootRadius = 0.3;
dustRadius = 4;

%% inputs into the I structure
I.elevation = p.Results.elevation;
I.mu0 = ucz;
I.SSA = ucssa;
I.muS = ucS;
I.dust = dust;
I.soot = soot;
I.bandPass = bandPass;
if exist('bands','var')
    I.bands = bands;
end

%% reflectance wavelengths use those at which refractive index is measured
[~,wavelength] = RefractiveIndex([],'ice','nm');
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

%N-D grids of the inputs for dirty snow
[V.elevation,V.mu0,V.SSA,V.muS,V.dust,V.soot] =...
    ndgrid(I.elevation,I.mu0,I.SSA,I.muS,I.dust,I.soot);

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
        subNames = {'SSA','muS','dust','soot'};
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
        lapArray = zeros(numel(subV.SSA),size(bandPass,1)).';
        loop_mu0 = I.mu0(kz);
        % this loop computes the clean-snow and dirty-snow values and then
        % the difference
        parfor kr=1:nLoop
            pv = W.Value;
            % clean diffuse reflectance, independent of illumination angle
            Pc = fwdPrescription('snow','cosZ',[],'SSA',...
                pv{1}.SSA(kr),'lookup',true,'sizeU','um',...
                'waveU','nm','wavelength',wavelength);
            [~,Rdif] =  SPIReS_fwd(Pc);
            %direct reflectance depends on local illumination angle
            Pc = fwdPrescription('snow','cosZ',pv{1}.muS(kr),'SSA',...
                pv{1}.SSA(kr),'lookup',true,'sizeU','um',...
                'waveU','nm','wavelength',wavelength);
            [~,Rb] =  SPIReS_fwd(Pc);
            reflTbl = table(wavelength,[Rb.reflectance Rdif.reflectance],...
                'VariableNames',{'wavelength','reflectance'});
            reflTbl.Properties.VariableUnits = {'nm',''};
            cleanBandRefl = bandPassReflectance(reflTbl,pscript,...
                'mu0',loop_mu0,'muS',pv{1}.muS(kr));
            % same calculation for dirty snow, diffuse
            Pd = fwdPrescription('snow','cosZ',[],'SSA',...
                pv{1}.SSA(kr),'lookup',true,'sizeU','um','waveU','nm',...
                'wavelength',wavelength,...
                'LAP',{'dust','soot'},'LAPfraction',[pv{1}.dust(kr) pv{1}.soot(kr)],...
                'LAPradius',[dustRadius sootRadius]);
            [~,Rdif] =  SPIReS_fwd(Pd);
            % direct
            Pd = fwdPrescription('snow','cosZ',pv{1}.muS(kr),'SSA',...
                pv{1}.SSA(kr),'lookup',true,'sizeU','um','waveU','nm',...
                'wavelength',wavelength,...
                'LAP',{'dust','soot'},'LAPfraction',[pv{1}.dust(kr) pv{1}.soot(kr)],...
                'LAPradius',[dustRadius sootRadius]);
            [~,Rb] =  SPIReS_fwd(Pd);
            reflTbl = table(wavelength,[Rb.reflectance Rdif.reflectance],...
                'VariableNames',{'wavelength','reflectance'});
            reflTbl.Properties.VariableUnits = {'nm',''};
            dirtyBandRefl = bandPassReflectance(reflTbl,pscript,...
                'mu0',loop_mu0,'muS',pv{1}.muS(kr));
            lapArray(:,kr) = (dirtyBandRefl-cleanBandRefl).'; %#ok<PFOUS>
        end
        % save as a checkpoint
        checkname = sprintf('checkDirty_%d_%d.mat',ka,kz);
        save(fullfile(folder,checkname),'ka','kz','lapArray','-v7.3')
        fprintf('loop elevation %d angle %d completed\n',ka,kz);
    end
end
files = dir(fullfile(folder,'checkDirty*.mat'));
end