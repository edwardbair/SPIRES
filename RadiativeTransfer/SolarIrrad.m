function [ Tvalues ] = SolarIrrad(varargin )
% [ Tbl ] = SolarIrrad(wavelength,Name/Value pairs)
% [Tbl = SolarIrrad(Name/Value pairs)
%
%Returns solar irradiance over the specified set of wavelengths, or if no
%wavelength argument, the wavelengths in the original data, for the
%specified atmospheric dataset or model. You can also specify a wavelength
%range and an integral value to scale to, an option that is mainly useful
%if you want to match a measured value from a broadband sensor.
% All inputs are optional, but wavelength must be first if specified. All
% other inputs are Name/Value pairs.
%   wavelength - scalar or vector (not matrix) of desired wavelengths (if
%   omitted, uses all wavelengths in original data, so you include the 2nd
%   output argument to retrieve the wavelengths)
% Other input, name-value pairs in any order (generally any unambiguous
%   abbreviation of 3 or more letters works)
%   'units' - wavelength units, default 'mum' (same as 'um'), other options
%       include any length unit (e.g. 'nm'), frequency (e.g. THz), or
%       wavenumber (cm^-1)
%   'atmosphere' - choices (case-insensitive, any 3+ letter abbreviation) are
%       all for clear atmospheres, but can be combined with a model for
%       transmittance through cloud, default 'SMARTS'
%       'ASTMG173' - ASTM standard for solar energy, sea level, U.S.
%           Standard atmosphere, adjusted to flat surface to match other
%           options
%       'SBDART' - Santa Barbara DISORT Atmospheric Radiative Transfer,
%           rural aerosol model, mid-latitude winter atmosphere, snow at
%           surface, 2000 m elevation, MODTRAN3 solar spectrum
%       'SMARTS' - Simple Model of the Atmospheric Radiative Transfer of
%           Sunshine, rural aerosol model, mid-latitude winter atmosphere,
%           snow at surface, 2000 m elevation, synthetic exoatmospheric
%           spectrum 1
%       (note that all these models go out only to 4 um, so values in the 4-5
%       um region are filled in with MODTRAN values and ATRAN transmission
%       estimates)
%   'solarZ' - solar zenith angle, degrees, default 48.19, so specify 0.0
%       to get the "solar constant"
%The following values are needed if you want to match the spectrum with an
%integral or even two or more integrals based on N multiple instruments
%   'waveBand' - matrix of size N x 2, output is scaled to match an
%       integrated value over this range (included to match broadband sensors),
%       no default, only used if 'integral' is specified
%   'integral' - scalar or vector of length N to match 'wavelengthRange',
%       no default so must be specified if 'wavelengthRange' is specified
%   'measType' - either scalar or cell vector of length N, specifying
%       'total', 'direct' (on flat surface), or 'diffuse'
%
% Output
%   Tvalues - table of values of solar irradiance under specified output conditions
%
% Reference: NED -- add the sources
% Gueymard, C. A. (2004), The sun's total and spectral irradiance for solar
% energy applications and solar radiation models, Solar Energy, 76, 423-453,
% doi: 10.1016/j.solener.2003.08.039.
%
% Values from 4-5 um filled in with MODTRAN values and the ATRAN
% transmission estimates:
% http://rredc.nrel.gov/solar/spectra/am0/modtran.html
% https://atran.sofia.usra.edu/cgi-bin/atran/atran.cgi

%% process inputs

persistent Fsmarts Fsbdart Fastm X already % lookup interpolation functions

p = inputParser;
% inputs and defaults
defaultSolarZ = 48.19;
mu0 = cosd(defaultSolarZ);
defaultUnits = 'mum';
addOptional(p,'lambda',[],@(x) isnumeric(x) && all(x(:)>0))
addParameter(p,validatestring('units',{'u','uni','units'}),defaultUnits,@ischar)
addParameter(p,validatestring('atmosphere',{'atm','atmos','atmosphere'}),...
    'smarts',@ischar)
addParameter(p,validatestring('solarz',{'sol','solarz'}),defaultSolarZ,...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<90)
% All of the following are needed if a spectrum is to match with an
% integral
addParameter(p,validatestring('waveband',{'waver','waverange','wavelengthrange','waveband'}),...
    [],@(x) isempty(x) || (isnumeric(x) && size(x,2)==2))
addParameter(p,validatestring('integral',{'int','integral'}),[],...
    @(x) isempty(x) || (isnumeric(x) && isvector(x) && all(x(:)>0)))
addParameter(p,validatestring('measurementtype',{'meas','type','meastype','measurementtype'}),...
    '',@(x) ischar(x) || iscell(x))

%parse input
parse(p,varargin{:});
wavelength = p.Results.lambda;
useTableWave = isempty(wavelength);
if ~useTableWave
    wavelength = unique(sort(wavelength(:)));
end
waveUnits = p.Results.units;
wavelengthRange = p.Results.waveband;
integratedValue = p.Results.integral;
measurementType = p.Results.measurementtype;
assert(~xor(isempty(wavelengthRange),isempty(integratedValue)) &&...
    ~xor(isempty(wavelengthRange),isempty(measurementType)),...
    'if ''wavelengthRange'' is specified, so must ''integral'' and ''measurementType''')
if ~isempty(wavelengthRange) && ~useTableWave
    warning('if you provide ''wavelengthRange'' you don''t need ''wavelength''')
end

atmosphere = validatestring(upper(p.Results.atmosphere),...
    {'SMARTS','SBDART','ASTMG173'});

% actual solar Z
solarZ = p.Results.solarz;

%% load values and build lookup functions
if isempty(already)
    already = true;
    X = load('solar_spectra.mat');
end
method = fittype('pchipinterp');
switch atmosphere
    % lookup functions are in order for the default condition:
    % 'toa','tau' (optical depth), 'xi' (factor for diffuse transmission)
    case 'SMARTS'
        if isempty(Fsmarts)
            % SMARTS values are for plane facing sun (apparently ... based on
            % values)
            x = X.smarts.WL;
            T = X.smarts;
            % beam radiation at top of atmosphere normal to sun
            Fsmarts{1} = fit(x,T.TOPDN,method);
            tau = -mu0*log(T.BOTDIR./T.TOPDN); % optical depth, independent of mu0
            Fsmarts{2} = fit(x,tau,method);
            Tdirect = exp(-tau/mu0);
            Tdif = T.BOTDIF./T.TOPDN;
            xi = mu0.*(1-Tdif-Tdirect)./tau;
            Fsmarts{3} = fit(x,xi,method); % xi factor for transmittance of diffuse
        end
        F = Fsmarts;
    case 'SBDART'
        if isempty(Fsbdart)
            % SBDART values are for a flat surface (apparently ... based on
            % values
            x = X.sbdart.WL;
            T = X.sbdart;
            % beam radiation at top of atmosphere normal to sun
            Fsbdart{1} = fit(x,T.TOPDN/mu0,method);
            tau = -mu0*log(T.BOTDIR./T.TOPDN); % optical depth
            Fsbdart{2} = fit(x,tau,method);
            Tdirect = exp(-tau/mu0);
            Tdif = (T.BOTDN-T.BOTDIR)./T.TOPDN;
            xi = mu0.*(1-Tdif-Tdirect)./tau;
            Fsbdart{3} = fit(x,xi,method); % xi factor for transmittance of diffuse
        end
        F = Fsbdart;
    case 'ASTMG173'
        if isempty(Fastm)
            % ASTMG173 values are for plane facing sun (apparently ... based on
            % values) so code is the same as for smarts
            x = X.ASTMG173.WL;
            T = X.ASTMG173;
            % beam radiation at top of atmosphere normal to sun
                      Fastm{1} = fit(x,T.TOPDN,method);
            tau = -mu0*log(T.BOTDIR./T.TOPDN); % optical depth, independent of mu0
            Fastm{2} = fit(x,tau,method);
            Tdirect = exp(-tau/mu0);
            Tdif = T.BOTDIF./(mu0*T.TOPDN);
            xi = mu0.*(1-Tdif-Tdirect)./tau;
            Fsmarts{3} = fit(x,xi,method); % xi factor for transmittance of diffuse
        end
        F = Fastm;
    otherwise
        error('''atmosphere'' ''%s'' not recognize') % should not reach
end

%% adjust transmittance based on solar Z
if useTableWave
    x = F{1}.p.breaks';
else
    x = convertLengthUnits(wavelength,waveUnits,defaultUnits);
end

tau = F{2}(x);
directTrans = exp(-tau/cosd(solarZ));
diffuseTrans = 1-directTrans-F{3}(x).*tau/cosd(solarZ);
solarTOA = F{1}(x)*cosd(solarZ);
% set NaNs outside range of table
t = x<min(F{1}.p.breaks) | x>max(F{1}.p.breaks);
if any(t)
    warning('%d wavelength values outside range of table, those values set to NaN',...
        nnz(t))
    solarTOA(t) = NaN;
    directTrans(t) = NaN;
    diffuseTrans(t) = NaN;
end

%% calculate values for the flat surface
surfaceDirect = solarTOA.*directTrans;
surfaceDiffuse = solarTOA.*diffuseTrans;
surfaceTotal = surfaceDirect+surfaceDiffuse;
%% match the integrals
if ~isempty(integratedValue)
    waveCell = cell(1,length(integratedValue));
    % convert wavelength range to default units and sort each row
    w = sort(convertLengthUnits(wavelengthRange,waveUnits,defaultUnits),2);
    assert(min(x)<=min(w(:)) && max(x)>=max(w(:)),...
        '''wavelengthRange'' outside limits of radiation data wavelengths')
    %process as cell vector
    if ischar(measurementType)
        measurementType = {measurementType};
    end
    thisIntegral = zeros(length(measurementType),1);
    for k=1:length(measurementType)
        switch measurementType{k}
            case 'total'
                integrandF = fit(F{1}.p.breaks',surfaceTotal,fittype);
                thisIntegral(k) = integrate(integrandF,w(k,2),w(k,1));
            case 'direct'
                integrandF = fit(F{1}.p.breaks',surfaceDirect,fittype);
                thisIntegral(k) = integrate(integrandF,max(w),min(w));
            case 'diffuse'
                integrandF = fit(F{1}.breaks,surfaceDiffuse,fittype);
                thisIntegral(k) = integrate(integrandF,max(w),min(w));
            otherwise
                error('''measurementType'' ''%s'' not recognized',measurementType{k})
        end
        % store just the wavelengths within the wavelength range
        k1 = find(F{1}.p.breaks<w(k,1),1,'last');
        k2 = find(F{1}.p.breaks>w(k,2),1,'first');
        if isempty(k1)
            k1 = 1;
        end
        if isempty(k2)
            k2 = length(F{1}.p.breaks);
        end
        waveCell{k} = F{1}.p.breaks(k1:k2)';
        x = (thisF.p.breaks(k1:k2))';
    end
    
    %more to add here to do the adjustment
end

% adjust units if input not same as table
if ~strcmpi(waveUnits,defaultUnits)
    wavelength = convertLengthUnits(x,defaultUnits,waveUints);
    solarTOA = convertLengthUnits(solarTOA,waveUnits,defaultUnits); % opposite direction!
    surfaceDirect = convertLengthUnits(surfaceDirect,waveUnits,defaultUnits);
    surfaceTotal = convertLengthUnits(surfaceTotal,waveUnits,defaultUnits);
    surfaceDiffuse = convertLengthUnits(surfaceDiffuse,waveUnits,defaultUnits);
else
    wavelength = x;
end

% variable names for table
tunits{1} = waveUnits;
tunits{2} = '';
tunits{3} = '';
tunits{4} = ['W m^-2 ' waveUnits '^-1'];
tunits{5} = tunits{4};
tunits{6} = tunits{4};
tunits{7} = tunits{4};
vname = {'WL','transDirect','transDiffuse','TOPDN','BOTDN','BOTDIR','BOTDIF'};
Tvalues = table(wavelength,directTrans,diffuseTrans,solarTOA,surfaceTotal,...
    surfaceDirect,surfaceDiffuse,'VariableNames',vname);
Tvalues.Properties.VariableUnits = tunits;

end