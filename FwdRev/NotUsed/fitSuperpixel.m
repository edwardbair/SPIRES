function [ostruct,stats,LP] = fitSuperpixel(unknowns,wavelength,spectralSP,topography,backReflectance,varargin)
% [ostruct,stats,LP] = fitSuperpixel(unknowns,wavelength,spectralSP,topography,backReflectance,name/value pairs)
%fit snow properties to individual superpixel
%
%Input
% unknowns - cell vector of unknowns, among 'radius', 'fSCA', 'wetness',
%   'waterEquivalent', 'dust', 'dustRadius', 'soot', 'sootRadius','muS'
% wavelength - vector in nm
% spectralSP - one row from spectralSP in consolidateSuperpixels
% topography - corresponding row
% backReflectance - corresponding reflectance of non-snow endmember
%Optional input, Name/Value pairs
% 'method' - 'lsq' (default), 'norm', or 'spec'
% 'atmos' - atmospheric transmission function
% 'diff' - diffuse fraction function
%
%Output
% ostruct - solved values of the unknowns
% stats - statistics of the fit
% LP - prescription of the fitted snow properties

p = inputParser;
addRequired(p,'unknowns',@iscell)
addRequired(p,'wavelength',@isnumeric)
addRequired(p,'spectralSP',@(x) isvector(x) && isnumeric(x))
addRequired(p,'topography',@(x) isvector(x) && isnumeric(x))
addRequired(p,'backReflectance',@(x) isvector(x) && isnumeric(x))
addParameter(p,validatestring('atmosphere',{'atm','atmos','atmosphere'}),...
    [],@(x) strcmpi(class(x),'cfit') ||...
    contains(class(x),'interpolant','IgnoreCase',true))
addParameter(p,'method','lsq',@ischar)
addParameter(p,validatestring('diffusefraction',{'diff','difff','diffusefraction'}),...
    [],@(x) strcmpi(class(x),'cfit') ||...
    contains(class(x),'interpolant','IgnoreCase',true))
parse(p,unknowns,wavelength,spectralSP,topography,backReflectance,varargin{:})

% check the optional parameters
if isempty(p.Results.atmosphere)
    warning('no function for atmospheric weight supplied')
end
if isempty(p.Results.diffusefraction)
    warning('no function for diffuse fraction supplied')
end

assert(~(any(contains(unknowns,'shade','IgnoreCase',true)) &&...
    any(contains(unknowns,'muS','IgnoreCase',true))),...
    'can''t have both ''shade'' and ''muS'' as unknowns')

% set functions for reflectance = f(wavelength) and R0 = f(wavelength)
Rfun = fit(wavelength(:),spectralSP(:),'pchip');
fR0 = fit(wavelength(:),backReflectance(:),'pchip');
% fit a 2nd reflectance function if shade included
if any(contains(unknowns,'shade','IgnoreCase',true))
    x = [min(wavelength) max(wavelength)];
    y = [0 0];
    F = fit(x',y','poly1');
    C = {fR0,F};
    fR0 = C;
    t = contains(unknowns,'shade','IgnoreCase',true);
    unknowns = unknowns(~t);
end

% nominal snow prescription
S = SnowCloudLimits;
doDust = any(contains(unknowns,'dust','IgnoreCase',true));
doSoot = any(contains(unknowns,'soot','IgnoreCase',true));
assert(~(doDust && doSoot),'generally cannot solve for both dust and soot')
% general prescription
P0 = fwdPrescription('snow','wavelength',wavelength,'waveu','nm',...
    'elevation',topography(1),'cosZ',topography(2),'muS',topography(3),...
    'viewF',topography(4),'radius',S.defaultSnowRadius,...
    'diffuseFraction',p.Results.diffusefraction,'R0',fR0);
if doDust
    P0.snow.cleanSnow = false;
    P0.snow.LAP = {'dust'};
    P0.snow.LAPfraction = mean(S.dust);
    P0.snow.LAPradius = S.defaultDustRadius;
elseif doSoot
    P0.snow.cleanSnow = false;
    P0.snow.LAP = {'soot'};
    P0.snow.LAPfraction = mean(S.soot);
    P0.snow.LAPradius = S.defaultSootRadius;
end

% output from inversion
[ostruct,stats,LP] = SPIReS_inv_spectral(Rfun,'snow',unknowns,P0,...
    'method',p.Results.method,'atmos',p.Results.atmosphere);
end