function [ file, totalLength ] = prepMieInterpolant(folder,medium,varargin)
% [ file, totalLength ] = prepMieInterpolant(folder,medium [,waterFraction])
%calculates and writes Mie scattering coefficients omega, g, Qext etc
%
%Input
%   folder - where results go (drive omitted if Azure file share)
%   medium - 'ice', 'water', 'dust', 'soot', or 'wetSnow'
%       (to be augmented for different kinds of dust)
%Optional input
%   waterFraction - typically 0.01 to 0.15, increments of 0.01, needed only
%       when medium is 'wetSnow'
%
%Output
%   file - .mat files where output written in table format with columns
%   wavelength (um), radius (um), omega, g, Qext

% number of items for effective radius in interpolant
Ice.nSize = 4*(91+31);
Water.nSize = round(Ice.nSize/2);
Dust.nSize = 4*31;
Soot.nSize = 4*31;

% parallel processing
useParallel = true;

%wavelengths to use based on the Warren-Brandt tabulation for ice and the
%Hale-Querry tabulation for water
waveUnits = 'um';
limits = [.2 200];
if strcmpi(medium,'water')
    [~,wavelength] = RefractiveIndex([],'water',waveUnits);
else
    [~,wavelength] = RefractiveIndex([],'ice',waveUnits);
end
firstN = find(wavelength<limits(1),1,'last');
lastN = find(wavelength>limits(2),1,'first');
wavelength = wavelength(firstN:lastN);

if nargin<2
    if isdeployed
        disp(['file = ' mfilename ' folder medium'])
    else
        disp(['file = ' mfilename '(folder, medium)'])
    end
    file = '';
    return
end

% input
p = inputParser;
S = SnowCloudLimits;
wLimits = S.wetSnow;
waterFracValidation = @(x) isnumeric(x) && isscalar(x) &&...
    x>=wLimits(1) && x<=wLimits(2);
addRequired(p,'folder',@ischar)
addRequired(p,'medium',@ischar)
addOptional(p,'waterfraction',[],waterFracValidation)
parse(p,folder,medium,varargin{:})

[outputFolder,diaryFolder] = identifyFolders(p.Results.folder);

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
        radius = linspace(sqrt(S.snowRadius(1)),sqrt(S.snowRadius(2)),Water.nSize);
    otherwise
        error('medium ''%s'' not recognized',medium)
end
if strcmpi(medium,'wetsnow')
    assert(~isempty(p.Results.waterfraction),...
        'if the medium is ''wetSnow'', the optional argument waterFraction is needed')
    wf = p.Results.waterfraction; % mass fraction
    if wf==0
        N = RefractiveIndex(wavelength,'ice',waveUnits);
    else
        iceDensity = 917;
        waterDensity = 1000;
        density = [iceDensity waterDensity];
        massFraction = [(1-wf) wf];
        volFraction = massFraction./density;
        sumVF = sum(volFraction,2);
        volFraction = volFraction./sumVF;
        % volume-weighted average of the dielectric functions
        N = sqrt(volFraction(1)*RefractiveIndex(wavelength,'ice',waveUnits).^2+...
            volFraction(2)*RefractiveIndex(wavelength,'water',waveUnits).^2);
    end
else
    N = RefractiveIndex(wavelength,medium,waveUnits);
end
radius = convertUnits(radius,S.unitsSize,waveUnits); % (should be same anyway)

% index of refraction of air - need to revise if we adapt this to water
% as the medium
nm = 1.000293;

% matrices of all possibilities
[W,R] = ndgrid(wavelength,radius.^2);
CRefIn = zeros(size(W),'like',N);
for r=1:size(W,1)
    CRefIn(r,:) = N(r);
end

totalLength = numel(W);
if isdeployed
    h = StartAzureDiary(mfilename,diaryFolder,['startIdx' num2str(1)]); %#ok<NASGU>
end

% output filename based on medium and first/last values]
interval = [num2str(1) '_' num2str(totalLength)];
if strcmpi(medium,'wetSnow')
    wfpart = num2str(round(wf*100));
    fileDesignator = ['_' medium '_' wfpart '.mat'];
else
    fileDesignator = ['_' medium '_' interval '.mat'];
end
file = fullfile(outputFolder,['Mie' fileDesignator]);

% these routines don't deal with real(CRefIn)<1 (need to fix this)
t = real(CRefIn)/nm<1;
if nnz(t)
    CRefIn(t) = complex(nm*(1+10*eps),imag(CRefIn(t)));
end
W = W(:);
R = R(:);
CRefIn = CRefIn(:);
x = 2*pi*R./W; % Mie parameter
useCAM = x>=20; % complex angular momentum, for large x

% compute Mie parameters and save to table
if strcmpi(medium,'wetSnow')
    variables = {'wavelength','radius','waterFraction','omega','g','Qext','Qabs','Qsca','Qpr','x','N'};
    waterF = repmat(wf,size(W));
else
    variables = {'wavelength','radius','omega','g','Qext','Qabs','Qsca','Qpr','x','N'};
end
if any(useCAM)
    M = mieapx(CRefIn(useCAM)/nm,x(useCAM),useParallel);
    if strcmpi(medium,'wetSnow')
        TCAM = table(W(useCAM),R(useCAM),waterF(useCAM),M.omega,M.g,M.Qext,...
            M.Qabs,M.Qsca,M.Qpr,...
            x(useCAM),CRefIn(useCAM)/nm,...
            'VariableNames',variables);
    else
        TCAM = table(W(useCAM),R(useCAM),M.omega,M.g,M.Qext,M.Qabs,M.Qsca,M.Qpr,...
            x(useCAM),CRefIn(useCAM)/nm,...
            'VariableNames',variables);
    end
end
if any(~useCAM)
    % convert radius and wavelength to meters
    lambda = convertUnits(W(~useCAM),waveUnits,'m');
    radius = convertUnits(R(~useCAM),waveUnits,'m');
    N = CRefIn(~useCAM);
    M = MieStuff(radius,lambda,N,nm,useParallel);
    if strcmpi(medium,'wetSnow')
        TBH = table(W(~useCAM),R(~useCAM),waterF(~useCAM),M.omega,M.g,M.Qext,...
            M.Qabs,M.Qsca,M.Qpr,...
            2*pi*radius./lambda,N/nm,...
            'VariableNames',variables);
    else
        TBH = table(W(~useCAM),R(~useCAM),M.omega,M.g,M.Qext,M.Qabs,M.Qsca,...
            M.Qpr,...
            2*pi*radius./lambda,N/nm,...
            'VariableNames',variables);
    end
end

% consolidate the tables
if any(useCAM) && any(~useCAM)
    T = [TBH; TCAM];
elseif any(useCAM)
    T = TCAM;
elseif any(~useCAM)
    T = TBH;
end
T = sortrows(T); %#ok<NASGU>
save(file,'T')
if isdeployed
    disp(['file ' file ' written'])
end
end