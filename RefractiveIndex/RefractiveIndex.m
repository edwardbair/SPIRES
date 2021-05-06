function [ N, varargout] = RefractiveIndex( wavelength, substance, varargin )
% [ N [, wavelengths] ] = RefractiveIndex( wavelength, substance [, units] )
% complex refractive index of ice or water, simple dust and soot
%
% INPUT
%   wave - wavelength(s), default in um (if empty, return values from
%   table)
%   substance - either 'ice', 'water', 'dust', or 'soot', or variants on
%   'iceW' - Warren 1984 (can add temperature dependence for wavelength>167 um)
%   'iceWB' - Warren & Brandt 2008;
%   'iceP' - Picard 2016 (returns same as 'iceWB' except for wavelengths the
%       Picard dataset covers, 0.32 to 0.60 um)
%   (can be character variable or categorical)
%
% OPTIONAL INPUT, in order
%   units for wavelength, default 'um', but 'nm', 'mum', 'mm', 'cm', 'm',
%       are supported
%   temperature in deg C only for 'iceW' and wavelengths > 167 um, default -5
%   (variability with temperature not implemented yet, so this variable
%   will be ignored)
%
% OUTPUT
%   N - complex refractive index (imaginary part positive)
% OPTIONAL OUTPUT
%   wavelengths - wavelengths from original database

% DATA SOURCES
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
% dust (just one dust, but more to add here)
% mostly from Skiles, S.M., Painter, T., & Okin, G.S. (2017). A method to
% retrieve the spectral complex refractive index and single scattering optical
% properties of dust deposited in mountain snow. Journal of Glaciology, 63,
% 133-147. https://doi.org/10.1017/jog.2016.126
% and from Zender  http://dust.ess.uci.edu/smn/smn_cgd_200610.pdf
%
% soot
% Bond, T.C., & Bergstrom, R.W. (2006). Light absorption by carbonaceous
% particles: An investigative review. Aerosol Science and Technology, 40, 27-67.
% https://doi.org/10.1080/02786820500421521

% method - for ice and water real parts against log wave, log imag parts
% against log wave
% for dust and soot, real and imaginary parts against wave

numarg = 2;
narginchk(numarg,4)
nargoutchk(0,2)

persistent W84r W84i WB08r WB08i P16r P16i HQ73r HQ73i

% make substances to check categorical variables
ice = categorical({'ice'});
water = categorical({'water'});
dust = categorical({'dust'});
soot = categorical({'soot'});
icew = categorical({'icew'});
icewb = categorical({'icewb'});
icep = categorical({'icep'});

% inputs
p = inputParser;
addRequired(p,'wavelength',@(x) isnumeric(x) && all(x(:)>0))
addRequired(p,'substance',@(x) ischar(x) || iscategorical(x))
addOptional(p,'units','um',@ischar)
addOptional(p,'iceTemp',-5,@(x) isnumeric(x) && isscalar(x))
parse(p,wavelength,substance,varargin{:})

if isempty(p.Results.units)
    units = 'um';
else
    units = p.Results.units;
end
iceTemp = p.Results.iceTemp; %#ok<NASGU>

% choose substance
if ischar(p.Results.substance)
    substance = categorical({lower(p.Results.substance)});
else
    substance = categorical({lower(char(p.Results.substance))});
end

% generic 'ice' as substance uses 'iceP'
if substance==ice
    substance = icep;
end

% 1st pass, calculate the interpolation functions
switch substance
    case water
        if isempty(HQ73i)
            load('WaterIndexRefraction.mat','WaterHaleQuerry73')
            [HQ73r,HQ73i] = fitValues(WaterHaleQuerry73);
        end
    case icep
        if isempty(P16i)
            m = matfile('IceIndexRefraction.mat');
            Picard16 = m.Picard16;
            WarrenBrandt08 = m.WarrenBrandt08;
            [Wr,Wi] = fitValues(WarrenBrandt08,'method','linear','extrap','linear','smooth',false);
            [~,P] = fitValues(Picard16,'method','linear','extrap','linear','smooth',false);
            compositeWave = unique([Wr.GridVectors{1}; P.GridVectors{1}]);
            t = compositeWave<=max(P.GridVectors{1});
            holdP = griddedInterpolant(compositeWave(t),P(compositeWave(t)));
            logImag = [max(Wi(compositeWave(t)),holdP(compositeWave(t)));...
                Wi(compositeWave(~t))];
            newTbl = table(exp(compositeWave),Wr(compositeWave),exp(logImag),...
                'VariableNames',{'wavelength','realpart','imagpart'});
            newTbl.Properties.VariableUnits = WarrenBrandt08.Properties.VariableUnits;
            [P16r,P16i] = fitValues(newTbl);
        end
    case icewb
        if isempty(WB08i)
            load('IceIndexRefraction.mat','WarrenBrandt08')
            [WB08r,WB08i] = fitValues(WarrenBrandt08);
        end
    case icew
        if isempty(W84i)
            m = matfile('IceIndexRefraction.mat');
            Warren84 = m.Warren84;
            if max(convertLengthUnits(p.Results.wavelength,units,Warren84.Properties.VariableUnits{1}))>...
                    max(Warren84.wavelength)
                Warren84microwave = m.Warren84microwave;
                % right now use -5C temperature, add temperature variability sometime
                holdW = Warren84microwave(:,{'wavelength','realpart5','imagpart5'});
                holdW.Properties.VariableNames = Warren84.Properties.VariableNames;
                Warren84 = [Warren84(1:end-1,:); holdW];
            end
            [W84r,W84i] = fitValues(Warren84);
        end
end

if ~isempty(wavelength)
    % wavelength in correct units if argument specified different than um
    wave = wavelength;
    wavelength = convertLengthUnits(wavelength,p.Results.units,'um');
    logWave = log(wavelength);
else
    switch substance
        case water
            [logWave,wave] = setWavelength(HQ73r,p.Results.units);
        case icep
            [logWave,wave] = setWavelength(P16r,p.Results.units);
        case icew
            [logWave,wave] = setWavelength(W84r,p.Results.units);
        case icewb
            [logWave,wave] = setWavelength(WB08r,p.Results.units);
        otherwise
            error('wavelength (first argument) cannot be empty unless substance is water or one of the ices')
    end
end

% all passes, do the interpolation
% if wavelength is a matrix, need to vectorize and then reshape results
if exist('logWave','var')
    waveSize = size(logWave);
    logWave = logWave(:);
    wavelength = exp(logWave);
else
    waveSize = size(wavelength);
    wavelength = wavelength(:);
end

% interpolation depending on substance
switch substance
    case icep
        Fr = P16r;
        Fi = P16i;
    case water
        Fr = HQ73r;
        Fi = HQ73i;
    case icew
        Fr = W84r;
        Fi = W84i;
    case icewb
        Fr = WB08r;
        Fi = WB08i;
    case dust
        [realpart,imagpart] = dustRefractive(wavelength);
    case soot
        [realpart,imagpart] = sootRefractive(wavelength);
    otherwise
        error(['substance ' substance ' not recognized'])
end
% calculate real and imaginary parts
if exist('Fr','var')
    realpart = Fr(logWave);
    imagpart = exp(Fi(logWave));
    if any(logWave<min(Fr.GridVectors{1})) ||...
            any(logWave>max(Fr.GridVectors{1}))
        warning('some wavelengths outside interpolation region, so ''nearest'' values returned')
    end
end

N = complex(realpart,abs(imagpart));
if ~isequal(size(N),waveSize)
    N = reshape(N,waveSize);
end

if nargout>1
    varargout{1} = wave;
end
end

function [n,k] = dustRefractive(lambda)
persistent dustAlready DI
if isempty(dustAlready)
    dustAlready = true;
    % data mostly from Skiles et al 2016, J Glac, doi 10.1017/jog.2016.126
    % and from Zender  http://dust.ess.uci.edu/smn/smn_cgd_200610.pdf
    % columns: wavelength (um), imag part, weight
    dustData = [0.35 0.0015 1; 0.40 0.0013 1; 0.45 0.0011 1; 0.50 0.0010 1;...
        0.55 0.0008 1; 0.60 0.0006 1; 0.65 0.0007 1; 0.70 0.0006 1;...
        0.75 0.0006 1; 0.80 0.0006 1; 0.90 0.0006 1; 1.00 0.0006 1;...
        1.10 0.0005 1; 1.20 0.0005 1; 1.30 0.0005 1; 1.40 0.0005 1;...
        1.50 0.0005 1; 0.30 0.0048 0.5; 0.35 0.0043 0.5; 0.40 0.0045 0.5;...
        0.45 0.0030 0.5; 0.50 0.0025 0.5; 0.55 0.0020 0.5; 0.60 0.0015 0.5;...
        0.65 0.0012 0.5; 0.70 0.0010 0.5; 0.75 0.0010 0.5; 0.80 0.0009 0.75;...
        0.90 0.0008 0.75; 1.00 0.0007 0.75; 1.10 0.0006 0.75; 1.80 0.0005 1;...
        1.90 0.0006 1; 2.00 0.0006 1; 2.10 0.0005 1; 2.20 0.0006 1;...
        2.30 0.0006 1; 2.40 0.0007 1; 2.50 0.0010 1];
    % smooth parameter from visual inspection of cftool results
    DI = fit(dustData(:,1),dustData(:,2),'smoothingspline',...
        'SmoothingParam',0.97,'Weights',dustData(:,3));
end
n = repmat(1.55,size(lambda));
if any(lambda<min(DI.p.breaks(1))) || any(lambda>max(DI.p.breaks(end)))
    warning('some wavelengths outside interpolation range, ''nearest'' returned')
end
lambda(lambda<DI.p.breaks(1)) = DI.p.breaks(1);
lambda(lambda>DI.p.breaks(end)) = DI.p.breaks(end);
k = DI(lambda);
end

function [n,k] = sootRefractive(lambda)
% value from Bond, T.C., and R.W. Bergstrom (2006), Light absorption by
% carbonaceous particles: An investigative review, Aerosol Science and
% Technology, 40, 27-67, doi: 10.1080/02786820500421521.
n = repmat(1.95,size(lambda));
k = repmat(0.79,size(lambda));
end

function [logW,W] = setWavelength(F,units)
logW = F.GridVectors{1};
W = convertLengthUnits(exp(logW),'um',units);
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
