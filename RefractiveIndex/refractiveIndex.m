function [ N, varargout] = refractiveIndex(wavelength,substance,waveunit  )
% [ N [,wavelengths]] = refractiveIndex( wavelength, substance, waveunit  )
% complex refractive index of specified substance
%
% INPUT
%   wavelength - (if empty, return values from table)
%   substance, case-insensitive, can be character or string
%       options for ice are:
%       'iceW' - Warren 1984 (can add temperature dependence for wavelength>167 um)
%       'iceWB' - Warren & Brandt 2008;
%       'iceP' - Picard 2016 (returns same as 'iceWB' except for wavelengths the
%           Picard dataset covers, 0.32 to 0.60 um)
%   'water'
%   Options for light-absorbing particles (first 5 characters enough) are:
%       'soot' - black carbon, density 1270 kg/m3
%       'coatsoot' - sulfur-coated black carbon
%       'brC' - brown carbon (Kirschstetter et al. 2004)
%       'coatbrC' - sulfur-coated brown carbon
%       'SaharaDust' (Balkanski et al. 2007)
%       'GreenlandDust' (Polashenski et al., 2015)
%       'SanJuanDust' (Skiles et al. 2017)
%   waveunit, no default to prevent errors, options are 'um','nm', 'mum',
%       'mm', 'cm', 'm',
%
% OUTPUT
%   N - complex refractive index (imaginary part positive)
% OPTIONAL OUTPUT
%   wavelengths - wavelengths from original database, in units specified by
%       waveunit, useful if input wavelengths = []

% DATA SOURCES
%
% ice
% Warren, S.G. (1984). Optical constants of ice from the ultraviolet to the
% microwave. Applied Optics, 23, 1206-1225. https://doi.org/10.1364/AO.23.001206

% Warren, S.G., & Brandt, R.E. (2008). Optical constants of ice from the
% ultraviolet to the microwave: A revised compilation. Journal of Geophysical
% Research: Atmospheres, 113, D14220. https://doi.org/10.1029/2007JD009744
% data from https://atmos.washington.edu/ice_optical_constants/
%
% Picard, G., Libois, Q., & Arnaud, L. (2016). Refinement of the ice absorption
% spectrum in the visible using radiance profile measurements in Antarctic snow.
% The Cryosphere, 10, 2655-2672. https://doi.org/10.5194/tc-10-2655-2016
% data from http://pp.ige-grenoble.fr/pageperso/picardgh/ice_absorption/
%
% water
% Hale, G.M., & Querry, M.R. (1973). Optical constants of water in the 200-nm
% to 200-µm wavelength region. Applied Optics, 12, 555-563.
% https://doi.org/10.1364/AO.12.000555
%
% black and brown carbon
% Bond, T.C., & Bergstrom, R.W. (2006). Light absorption by carbonaceous
% particles: An investigative review. Aerosol Science and Technology, 40, 27-67.
% https://doi.org/10.1080/02786820500421521
%
% Kirchstetter, T. W., T. Novakov, and P. V. Hobbs (2004), Evidence
% that the spectral dependence of light absorption by aerosols is
% affected by organic carbon, J. Geophys. Res., 109, D21208,
% https://doi.org/10.1029/2004JD004999
%
% Dust
% Balkanski, Y., M. Schulz, T. Claquin, and S. Guibert (2007),
% Reevaluation of mineral aerosol radiative forcings suggests a better
% agreement with satellite and AERONET data, Atmos. Chem. Phys., 7
% (1), 81-95, https://doi.org/10.5194/acp-7-81-2007.
%
% Polashenski, C. M., J. E. Dibb, M. G. Flanner, J. Y. Chen,
% Z. R. Courville, A. M. Lai, J. J. Schauer, M. M. Shafer, and
% M. Bergin (2015), Neither dust nor black carbon causing apparent
% albedo decline in Greenland's dry snow zone: Implications for MODIS
% C5 surface reflectance, Geophys. Res. Lett., 42, 9319-9327,
% https://doi.org/10.1002/2015GL065912.
%
% Skiles, S. M., Painter, T., and Okin, G. (2017), A method to
% retrieve the spectral complex refractive index and single scattering
% optical properties of dust deposited in mountain snow, Journal of
% Glaciology. 63, 133-147, https://doi.org/10.1017/jog.2016.126.

% method - real parts against log wave, log imag parts against log wave

numarg = 3;
narginchk(numarg,4)
nargoutchk(0,2)

persistent W84 WB08 P16 HQ73 BC CBC BrC CBrC Sahara Greenland SanJuan alreadyWarn

tableUnit = 'm';

% inputs
p = inputParser;
addRequired(p,'wavelength',@(x) isnumeric(x) && all(x(:)>0))
addRequired(p,'substance',@(x) ischar(x) || isstring(x))
addRequired(p,'waveunit',@(x) ischar(x) || isstring(x))
parse(p,wavelength,substance,waveunit);
% iceTemp = p.Results.iceTemp; %#ok<NASGU>

% change substance to a unique categorical variable
substance = assignSubstance(p.Results.substance);
FileName = fname(substance);

% make categorical variables to check against substance
water = categorical({'water'});
iceW = categorical({'iceW'});
iceWB = categorical({'iceWB'});
iceP = categorical({'iceP'});
BlackCarbon = categorical({'BlackCarbon'});
CoatedBlackCarbon = categorical({'CoatedBlackCarbon'});
BrownCarbon = categorical({'BrownCarbon'});
CoatedBrownCarbon = categorical({'CoatedBrownCarbon'});
SaharaDust = categorical({'SaharaDust'});
GreenlandDust = categorical({'GreenlandDust'});
SanJuanDust = categorical({'SanJuanDust'});

% 1st pass for any substance, calculate the interpolation functions
switch substance
    case water
        if isempty(HQ73)
            m = matfile(FileName);
            [HQ73{1},HQ73{2}] = fitValues(m.HaleQuerry73);
        end
    case iceP
        if isempty(P16)
            m = matfile(FileName);
            WarrenBrandt08 = m.WarrenBrandt08;
            [Wr,Wi] = fitValues(m.WarrenBrandt08,'method','linear','extrap','nearest','smooth',false);
            [~,P] = fitValues(m.Picard16,'method','linear','extrap','linear','smooth',false);
            compositeWave = unique([Wr.GridVectors{1}; P.GridVectors{1}]);
            t = compositeWave<=max(P.GridVectors{1});
            holdP = griddedInterpolant(compositeWave(t),P(compositeWave(t)));
            logImag = [max(Wi(compositeWave(t)),holdP(compositeWave(t)));...
                Wi(compositeWave(~t))];
            newTbl = table(exp(compositeWave),Wr(compositeWave),exp(logImag),...
                'VariableNames',{'wavelength','realpart','imagpart'});
            newTbl.Properties.VariableUnits = WarrenBrandt08.Properties.VariableUnits;
            [P16{1},P16{2}] = fitValues(newTbl);
        end
    case iceWB
        if isempty(WB08)
            m = matfile(FileName);
            [WB08{1},WB08{2}] = fitValues(m.WarrenBrandt08);
        end
    case iceW
        if isempty(W84)
            m = matfile(FileName);
            Warren84 = m.Warren84;
            if max(convertLengthUnits(p.Results.wavelength,waveunit,Warren84.Properties.VariableUnits{1}))>...
                    max(Warren84.wavelength)
                Warren84microwave = m.Warren84microwave;
                % right now use -5C temperature, add temperature variability sometime
                holdW = Warren84microwave(:,{'wavelength','realpart5','imagpart5'});
                holdW.Properties.VariableNames = Warren84.Properties.VariableNames;
                Warren84 = [Warren84(1:end-1,:); holdW];
            end
            [W84{1},W84{2}] = fitValues(Warren84);
        end
    case BlackCarbon
        if isempty(BC)
            m = matfile(FileName);
            [BC{1},BC{2}] = fitValues(m.Soot);
        end
    case CoatedBlackCarbon
        if isempty(CBC)
            m = matfile(FileName);
            [CBC{1},CBC{2}] = fitValues(m.CoatedSoot);
        end
    case BrownCarbon
        if isempty(BrC)
            m = matfile(FileName);
            [BrC{1},BrC{2}] = fitValues(m.BrownCarbon);
        end
    case CoatedBrownCarbon
        if isempty(CBrC)
            m = matfile(FileName);
            [CBrC{1},CBrC{2}] = fitValues(m.CoatedBrownCarbon);
        end
    case SaharaDust
        if isempty(Sahara)
            m = matfile(FileName);
            [Sahara{1},Sahara{2}] = fitValues(m.SaharaDust);
        end
    case GreenlandDust
        if isempty(Greenland)
            m = matfile(FileName);
            [Greenland{1},Greenland{2}] = fitValues(m.GreenlandDust);
        end
    case SanJuanDust
        if isempty(SanJuan)
            m = matfile(FileName);
            [SanJuan{1},SanJuan{2}] = fitValues(m.SanJuanDust);
        end
end

% all passes, interpolation depending on substance
switch substance
    case iceP
        Fr = P16{1};
        Fi = P16{2};
    case water
        Fr = HQ73{1};
        Fi = HQ73{2};
    case iceW
        Fr = W84{1};
        Fi = W84{2};
    case iceWB
        Fr = WB08{1};
        Fi = WB08{2};
    case SaharaDust
        Fr = Sahara{1};
        Fi = Sahara{2};
    case GreenlandDust
        Fr = Greenland{1};
        Fi = Greenland{2};
    case SanJuanDust
        Fr = SanJuan{1};
        Fi = SanJuan{2};
    case BlackCarbon
        Fr = BC{1};
        Fi = BC{2};
    case CoatedBlackCarbon
        Fr = CBC{1};
        Fi = CBC{2};
    case BrownCarbon
        Fr = BrC{1};
        Fi = BrC{2};
    case CoatedBrownCarbon
        Fr = CBrC{1};
        Fi = CBrC{2};
    otherwise
        error(['substance ' char(substance) ' not recognized'])
end

if ~isempty(wavelength)
    % convert wavelength to correspond with table
    wavelength = convertLengthUnits(wavelength,p.Results.waveunit,tableUnit);
    logWave = log(wavelength);
else
    logWave = Fr.GridVectors{1};
end

% all passes, do the interpolation
% if wavelength is a matrix, need to vectorize and then reshape results
waveSize = size(logWave);
logWave = logWave(:);

% calculate real and imaginary parts
realpart = Fr(logWave);
imagpart = exp(Fi(logWave));
if any(logWave<min(Fr.GridVectors{1})) ||...
        any(logWave>max(Fr.GridVectors{1}))
    if ~alreadyWarn
        warning('some wavelengths outside interpolation region, so ''nearest'' values returned')
        alreadyWarn = true;
    end
end

N = complex(realpart,abs(imagpart));
if ~isequal(size(N),waveSize)
    N = reshape(N,waveSize);
end

if nargout>1
    wave = convertLengthUnits(exp(logWave),tableUnit,p.Results.waveunit);
    varargout{1} = reshape(wave,size(N));
end
end

function [realArg,imagArg] = fitValues(Tbl,varargin)
%fit to real and imagingary parts of refractive index

p = inputParser;
addRequired(p,'Tbl',@istable)
addParameter(p,'method','makima',@ischar)
addParameter(p,'extrap','nearest',@ischar)
addParameter(p,'smooth',true,@(x) islogical(x) || isnumeric(x))
parse(p,Tbl,varargin{:})
Tbl = p.Results.Tbl;
method = p.Results.method;
extrap = p.Results.extrap;
doSmooth = logical(p.Results.smooth);

logWave = log(Tbl.wavelength);
if any(contains(Tbl.Properties.VariableNames,'realpart'))
    if doSmooth
        F = fit(logWave,Tbl.realpart,'smoothingspline');
        realArg = griddedInterpolant(logWave,F(logWave),method,extrap);
    else
        realArg = griddedInterpolant(logWave,Tbl.realpart,method,extrap);
    end
else
    realArg = [];
end

if any(contains(Tbl.Properties.VariableNames,'imagpart'))
    if doSmooth
        F = fit(logWave,log(Tbl.imagpart),'smoothingspline');
        imagArg = griddedInterpolant(logWave,F(logWave),method,extrap);
    else
        imagArg = griddedInterpolant(logWave,log(Tbl.imagpart),method,extrap);
    end
else
    imagArg = [];
end
end

function material=assignSubstance(substance)
% assign substance to a unique categorical variable

% make categorical variables to check against substance
water = categorical({'water'});
iceW = categorical({'iceW'});
iceWB = categorical({'iceWB'});
iceP = categorical({'iceP'});
BlackCarbon = categorical({'BlackCarbon'});
CoatedBlackCarbon = categorical({'CoatedBlackCarbon'});
BrownCarbon = categorical({'BrownCarbon'});
CoatedBrownCarbon = categorical({'CoatedBrownCarbon'});
SaharaDust = categorical({'SaharaDust'});
GreenlandDust = categorical({'GreenlandDust'});
SanJuanDust = categorical({'SanJuanDust'});

% generic 'ice' as substance uses 'iceP'
if strcmpi(substance,'ice')
    material = iceP;
elseif strcmpi(substance,'icep')
    material = iceP;
elseif strcmpi(substance,'water')
    material = water;
elseif strcmpi(substance,'icew')
    material = iceW;
elseif strcmpi(substance,'icewb')
    material = iceWB;
elseif strcmpi(substance,'soot') || strncmpi(substance,'black',5)
    material = BlackCarbon;
elseif strncmpi(substance,'coats',5)
    material = CoatedBlackCarbon;
elseif strcmpi(substance,'brc')
    material = BrownCarbon;
elseif strncmpi(substance,'coatb',5)
    material = CoatedBrownCarbon;
elseif strncmpi(substance,'sahar',5)
    material = SaharaDust;
elseif strncmpi(substance,'green',5)
    material = GreenlandDust;
elseif strncmpi(substance,'sanju',5)
    material = SanJuanDust;
else
    error('substance %s not recognized',substance)
end
end

function fn=fname(substance)
% file name that has the values for the substance
if contains(char(substance),'water','IgnoreCase',true)
    fn = 'WaterIndexRefraction.mat';
elseif contains(char(substance),'ice','IgnoreCase',true)
    fn = 'IceIndexRefraction.mat';
elseif contains(char(substance),'dust','IgnoreCase',true)
    fn = 'DustIndexRefraction.mat';
elseif contains(char(substance),'carbon','IgnoreCase',true)
    fn = 'CarbonIndexRefraction.mat';
else
    error('can''t match file name for substance ''%s''',char(substance))
end
end
