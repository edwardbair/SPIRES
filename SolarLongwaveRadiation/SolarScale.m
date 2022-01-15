function [ Tbl ] = SolarScale(varargin )
% [ Tbl ] = SolarScale(Name, Value)
%Returns scaled solar radiation such that the integral from 0.28 to 4.0 um
%is 1.0.
%Therefore, to get spectral distribution of TOA or surface irradiance,
%multiply the result from SolarScale by the integrated value (i.e., exoatmospheric
%solar radiation, perhaps adjusted for radius vector and solar zenith angle
%if you want TOA spectrum)
%Values are based on the SMARTS298 model for a mid-latitude winter (other
%atmospheres to be added)
%
%Required software: the Shape Language Modeling Toolbox from the MATLAB
%File Exchange
%
%Required input
%Input, name-value pairs, can be scalars, vectors, or matrices,
%   but any that are non-scalar must be same size (any abbreviation of 3 or
%   more letters works)
%   You can use ndgrid to create inputs of the same size
% 'location' - 'TOA' or 'surface' (default)
% 'wavelength' - desired wavelengths (if omitted, uses all wavelengths
%   in original data)
% 'smooth' - if true, smooths spectrum using SLM toolbox
% 'units' - units for wavelength, no default to prevent errors
% 'cosZ' - cosine of solar zenith angle, default cosd(48.16)
% 'elevation' - in meters, default 3000
%
%Output
% Returns table of wavelengths with scaled values
%   of solar radiation, so that integral from 0.28 to 4 um is 1.0 (so one can
%   multiply by a total to get a spectrum), and diffuse fraction at each wavelength
%
%
%References
% Gueymard, C. A. (2018). Revised composite extraterrestrial spectrum based
%   on recent solar irradiance observations. Solar Energy, 169, 434-440.
%   https://doi.org/10.1016/j.solener.2018.04.067
% Gueymard, C. A. (2019). The SMARTS spectral irradiance model after 25 years:
%   New developments and validation of reference spectra. Solar Energy, 187,
%   233-253. https://doi.org/10.1016/j.solener.2019.05.048


% Future modification to fill in values for 4-5 um with MODTRAN values and
% the ATRAN transmission estimates:
% http://rredc.nrel.gov/solar/spectra/am0/modtran.html
% https://atran.sofia.usra.edu/cgi-bin/atran/atran.cgi

persistent FG FD FTOA lookup_waveU lookup_elevU

p = inputParser;
% inputs and defaults
addParameter(p,validatestring('units',{'unit','units'}),'',...
    @(x) (ischar(x) || isstring(x) || iscategorical(x)))
addParameter(p,validatestring('location',{'loc','loca','location'}),...
    'surface',@(x) ischar(x) || isstring(x) || iscategorical(x))
addParameter(p,'cosz',cosd(48.19),@(x) isnumeric(x) && all(x(:)>=0))
addParameter(p,validatestring('elevation',{'elev','elevation'}),3000,@isnumeric)
addParameter(p,validatestring('wavelength',{'wave','wvl','wavelength'}),...
    [],@(x) isnumeric(x) && all(x(:)>0))
addParameter(p,'smooth',false,@(x) isscalar(x) && (islogical(x) || isnumeric(x)))
parse(p,varargin{:});

% check units
assert(~isempty(p.Results.units),...
    '''units'' must be specified, typically ''nm'' or ''um'' but any length unit works')
waveUnits = char(p.Results.units);

% set options
doSurface = strcmpi(char(p.Results.location),'surface');
doTOA = strcmpi(char(p.Results.location),'toa');
assert(doSurface || doTOA,...
    '''location'' must be either ''surface'' or ''TOA''')
doSmooth = logical(p.Results.smooth);

% first pass, get interpolants
if isempty(FG) || isempty(FD) || isempty(FTOA)
    [FG,FD,FTOA,lookup_waveU,lookup_elevU] = getSolarScaleInterpolants('mlw');
end

%check sizes, convert elevation to km
[cosZ,elevation] = checkSizes(p.Results.cosz,...
    convertLengthUnits(p.Results.elevation,'m',lookup_elevU));

if doSmooth || isempty(p.Results.wavelength)
    % get values from interpolant, initially use all wavelengths because if we
    % want to smooth later, we use a different set
    [wavelength,cosZ,elevation] = ndgrid(FG.GridVectors{1},cosZ,elevation);
else
    % convert wavelength to nm
    [wavelength,cosZ,elevation] =...
        ndgrid(convertLengthUnits(p.Results.wavelength,p.Results.units,lookup_waveU),...
        cosZ,elevation);
end

if doSurface
    Global = FG(wavelength,cosZ,elevation);
    DiffuseFraction = FD(wavelength,cosZ,elevation);
    Tbl = table(wavelength(:),cosZ(:),1000*elevation(:),Global(:),...
        DiffuseFraction(:),'VariableNames',...
        {'wavelength','cosZ','elevation','Global','DiffuseFraction'});
elseif doTOA
    TOA = FTOA(wavelength);
    Tbl = table(wavelength(:),TOA(:),'VariableNames',...
        {'wavelength','TOA'});
else
    error('should not reach this statement, check code')
end

% convert units
if ~strcmpi(waveUnits,'nm')
    Tbl.wavelength = convertLengthUnits(Tbl.wavelength,lookup_waveU,waveUnits);
    % arguments for conversion reversed when values are per wavelength
    if doSurface
        Tbl.Global = convertLengthUnits(Tbl.Global,waveUnits,lookup_waveU);
        Tbl.Properties.VariableUnits = {waveUnits,'','m','',''};
    else
        Tbl.TOA = convertLengthUnits(Tbl.TOA,waveUnits,lookup_waveU);
        Tbl.Properties.VariableUnits = {waveUnits,''};
    end
else
    if doSurface
        Tbl.Properties.VariableUnits = {waveUnits,'','m','',''};
    else
        Tbl.Properties.VariableUnits = {waveUnits,''};
    end
end

% revise table to smooth with SLM (this assumes that input wavelengths are
% in the right units specified by 'units')
if ~isempty(p.Results.wavelength) && doSmooth
    % knot at every 5th wavelength input wavelength and at the ends
    P = slmset('knots',round(length(p.Results.wavelength)/5));
    if doTOA
        P = slmset(P,'leftvalue',Tbl.TOA(1),'rightvalue',Tbl.TOA(end));
        slm = slmengine(Tbl.wavelength,Tbl.TOA,P);
        newTbl = table(slm.knots,slmeval(slm.knots,slm),'VariableNames',...
            {'wavelength','TOA'});
        newTbl.Properties.VariableUnits = {waveUnits,''};
    else
        [ucos,uelev] = ndgrid(unique(Tbl.cosZ),unique(Tbl.elevation));
        ucos = ucos(:);
        uelev = uelev(:);
        newTbl = table;
        for k=1:length(ucos)
            t = Tbl.cosZ==ucos(k) & Tbl.elevation==uelev(k);
            n1 = find(t,1,'first');
            n2 = find(t,1,'last');
            P = slmset('leftvalue',Tbl.Global(n1),'rightvalue',Tbl.Global(n2));
            slm1 = slmengine(Tbl.wavelength(t),Tbl.Global(t),P);
            P = slmset(P,'leftvalue',Tbl.DiffuseFraction(n1),'rightvalue',...
                Tbl.DiffuseFraction(n2));
            slm2 = slmengine(Tbl.wavelength(t),Tbl.DiffuseFraction(t),P);
            thisTbl = table(slm1.knots,repmat(ucos(k),length(slm1.knots),1),...
                repmat(uelev(k),length(slm1.knots),1),slmeval(slm1.knots,slm1),...
                slmeval(slm2.knots,slm2),'VariableNames',...
                {'wavelength','cosZ','elevation','Global','DiffuseFraction'});
            newTbl = [newTbl;thisTbl]; %#ok<AGROW>
        end
        newTbl.Properties.VariableUnits = {waveUnits,'','m','',''};
    end
    Tbl = newTbl;
end

end