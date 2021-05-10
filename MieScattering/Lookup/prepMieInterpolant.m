function [ C, M  ] = prepMieInterpolant(medium,varargin)
% [ Cstruct, Mstruct ] = prepMieInterpolant(medium [,waterFraction])
%calculates Mie scattering coefficients omega, g, Qext etc to be used in
%creating lookup functions
%
%Input
%   medium - 'ice', 'water', 'dust', 'soot', or 'wetSnow'
%       (to be augmented for different kinds of dust)
%Optional input
%   waterFraction - typically 0.01 to 0.15, increments of 0.01, needed only
%       when medium is 'wetSnow'
%
%Output
%   Cstruct - structure with inputs as N-dimensional grids
%   Mstruct - Mie structure with outputs as N-D grids

% number of items for effective radius in interpolant
Ice.nSize = 4*(91+31);
Water.nSize = 8*31;
Dust.nSize = 4*31;
Soot.nSize = 4*31;
WetSnow.nSize = 31;

%wavelengths to use based on the Picard and Warren-Brandt tabulations for ice
%and the Hale-Querry tabulation for water
waveUnits = 'um';
limits = [.28 20];
if strcmpi(medium,'water')
    [~,wavelength] = RefractiveIndex([],'water',waveUnits);
else
    [~,wavelength] = RefractiveIndex([],'ice',waveUnits);
end
firstN = find(wavelength<limits(1),1,'last');
lastN = find(wavelength>limits(2),1,'first');
wavelength = wavelength(firstN:lastN);

% input
p = inputParser;
S = SnowCloudLimits;
wLimits = S.wetSnow;
waterFracValidation = @(x) isnumeric(x) &&...
    all(x>=wLimits(1)) && all(x<=wLimits(2));
addRequired(p,'medium',@ischar)
addOptional(p,'waterFraction',[],waterFracValidation)
parse(p,medium,varargin{:})

% size intervals based on medium
switch lower(medium)
    case 'ice'
        radius = linspace(sqrt(S.iceCloudRadius(1)),sqrt(S.snowRadius(2)),Ice.nSize);
    case 'water'
        radius = linspace(sqrt(S.waterCloudRadius(1)),sqrt(S.waterInSnowRadius(2)),Water.nSize);
    case 'dust'
        radius = linspace(sqrt(S.dustRadius(1)),sqrt(S.dustRadius(2)),Dust.nSize);
    case 'soot'
        radius = linspace(sqrt(S.sootRadius(1)),sqrt(S.sootRadius(2)),Soot.nSize);
    case 'wetsnow'
        radius = linspace(sqrt(S.snowRadius(1)),sqrt(S.snowRadius(2)),WetSnow.nSize);
    otherwise
        error('medium ''%s'' not recognized',medium)
end

% matrices of all wavelength-radius combinations
[W,R] = ndgrid(wavelength,radius.^2);

% refractive index
if strcmpi(medium,'wetsnow')
    assert(~isempty(p.Results.waterFraction),...
        'if the medium is ''wetSnow'', the optional argument waterFraction is needed')
    wf = p.Results.waterFraction(:); % mass fraction
    iceDensity = 917;
    waterDensity = 1000;
    density = repmat([iceDensity waterDensity],length(wf),1);
    massFraction = [(1-wf) wf];
    volFraction = massFraction./density;
    sumVF = sum(volFraction,2);
    volFraction = volFraction./sumVF;
else
    CRefIn = RefractiveIndex(W,medium,waveUnits);
end
R = convertLengthUnits(R,S.unitsSize,waveUnits); % (should be same anyway)

% index of refraction of air - need to revise if we adapt this to water
% as the medium
nm = 1.000293;
x = 2*pi*R./W; % Mie parameter

% Mie values
if ~strcmpi(medium,'wetsnow')
    CRefIn = complex(real(CRefIn),abs(imag(CRefIn))); % make sure imag part positive
    M = MieSphere(R,'um',W,'um','refindex',CRefIn, 'mediumRI',nm);
    C.wavelength = W;
    C.radius = R;
    C.x = x;
    C.CRefIn = CRefIn;
    C.mediumRI = nm;
else
    [waterF,wW,wR] = ndgrid(wf,wavelength,radius.^2);
    mvar = {'omega','g','Qext'};
    for v=1:length(mvar)
        M.(mvar{v}) = zeros(size(waterF));
    end
    % volume-weighted average of the dielectric functions
    for k=1:length(wf)
        t = waterF==wf(k);
        CRefIn = sqrt(volFraction(k,1)*...
            RefractiveIndex(W,'ice',waveUnits).^2+...
            volFraction(k,2)*...
            RefractiveIndex(W,'water',waveUnits).^2);
        CRefIn = complex(real(CRefIn),abs(imag(CRefIn)));
        holdM = MieSphere(R,'um',W,'um','refindex',CRefIn,'mediumRI',nm);
        for v=1:length(mvar)
            M.(mvar{v})(t) = holdM.(mvar{v});
        end
    end
    C.wavelength = wW;
    C.radius = wR;
    C.waterF = waterF;
    C.mediumRI = nm;
end
end