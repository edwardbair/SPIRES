function [Ftrans,Fdf] = buildAtmosWeightF(wavelength,varargin)
% [Ftrans,Fdf] = buildAtmosWeightF(wavelength, name/value pairs)
%
% build functions for spectral atmospheric weighting and fraction of
% incoming irradiance that is diffuse, used to weight spectroscopic data
% for subsequent analyses
%
% requires SMARTS atmospheric model and SLM shape language modeling
% toolbox, but functions once created do not need them
%
% Input
%   wavelength - vector of wavelengths for the output functions
% Optional input, name/value pairs in either order
%   'atmos' - 3-letter atmospheric code, case-insensitive, default 'MLW'
%   'waveUnit' - default 'nm'
%
% Output, 2 interpolants
%   Ftrans - atmospheric transmittance, 0-1, (direct + diffuse/2)^2 as a function
%       of wavelength (squared because of transmission in both directions)
%   Fdf - diffuse fraction of irradiance as a function of wavelength
%
%References
% Gueymard, C. A. (2018). Revised composite extraterrestrial spectrum based
%   on recent solar irradiance observations. Solar Energy, 169, 434-440.
%   https://doi.org/10.1016/j.solener.2018.04.067
% Gueymard, C. A. (2019). The SMARTS spectral irradiance model after 25 years:
%   New developments and validation of reference spectra. Solar Energy, 187,
%   233-253. https://doi.org/10.1016/j.solener.2019.05.048

%% inputs
p = inputParser;
addRequired(p,'wavelength',@(x) isnumeric(x) && isvector(x) && all(x>0))
addParameter(p,'atmos','MLW',@ischar)
addParameter(p,validatestring('waveunit',{'wu','wave','waveu','waveunit'}),...
    'nm',@ischar)
parse(p,wavelength,varargin{:})
wavelength = unique(p.Results.wavelength(:)); % just in case

%% fixed parameter ranges
altit = (1000:1000:4000);
cosZ = linspace(.15,1,6);
[w,Z,mu0] = ndgrid(wavelength,altit,cosZ);
Trans = zeros(size(w));
DiffuseFraction = zeros(size(w));

%% SMARTS values
SetSMARTSversion(298);
for n=1:length(altit)
    for m=1:length(cosZ)
        P = defaultSMARTSinput(p.Results.atmos,'altit',altit(n),'cosZ',cosZ(m));
        S = SMARTSMain(P);
        % raw SMARTS output at native wavelengths, nm
        thisT = S.spectralTbl(:,{'waveL','HorzDiffuse','HorzGlobal','DirectTrans','TransDiffuse'});
        if ~strcmpi(p.Results.waveunit,'nm')
            thisT.waveL = convertLengthUnits(thisT.waveL,'nm',p.Results.waveunit);
        end
        % find the knots values for all passes through loops
        if n==1 && m==1
            k1 = find(thisT.waveL<wavelength(1),1,'last');
            k2 = find(thisT.waveL>wavelength(end),1,'first');
            knots = [thisT.waveL(k1); wavelength; thisT.waveL(k2)];
        end
        thisT = thisT(k1:k2,:);
        slmT = slmengine(thisT.waveL,(thisT.DirectTrans+thisT.TransDiffuse/2).^2,'knots',knots);
        slmD = slmengine(thisT.waveL,thisT.HorzDiffuse./thisT.HorzGlobal,'knots',knots);
        % output values
        k = Z==altit(n) & mu0==cosZ(m);
        Trans(k) = slmeval(w(k),slmT);
        DiffuseFraction(k) = slmeval(w(k),slmD);
    end
end
% rescale to range .001-.999 (not 0-1 so avoid under/overshoot
Trans = rescale(Trans,.001,.999);
DiffuseFraction = rescale(DiffuseFraction,.001,.999);
Ftrans = griddedInterpolant(w,Z,mu0,Trans,'linear','nearest');
Fdf = griddedInterpolant(w,Z,mu0,DiffuseFraction,'linear','nearest');

end