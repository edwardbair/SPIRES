function [Refl, Trans, Beam] = twostream( mu0, omega, g, method, varargin )
% [Refl, Trans, Beam] = twostream( mu0, omega, g, method [,Name/Value] )
%twostream: Two-stream solution of radiative transfer for single layer
%
% Input (omega & g must be same size, mu0 can be scalar or same size)
%   mu0 - cosine of illumination angle (enter empty value to calculate
%       reflectance to diffuse illumination)
%   omega - single-scattering albedo
%   g - asymmetry factor
%   method - either 'delta-Eddington' or 'Meador-Weaver hybrid'
%       [any 3- or more-letter abbreviation works]
% Optional input - name value pairs
%   'tau' - optical depth of the layer, must be scalar or same size as omega & g,
%       default is Inf
%   'R0', then number - reflectance of substrate, set to zero if omitted
%       (scalar or same size as omega & g)
%   'airmass', then number - compensated for Earth curvature, set to 1/mu0
%       if omitted, otherwise must be same size as mu0, see kasten.m for
%       calculations in Earth's atmosphere but omit for modeling snow
%       albedo
%
% Output
%   Refl - Reflectance of layer, including substrate
%   Trans - Diffuse ransmittance of layer
%   Beam - Direct transmittance of layer
%
% sources:
% Meador, W. E., & Weaver, W. R. (1980). Two-stream approximations to radiative
% transfer in planetary atmospheres – A unified description of existing methods
% and a new improvement. Journal of the Atmospheric Sciences, 37, 630-643.
% doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
% Wiscombe, W. J., & Warren, S. G. (1980). A model for the spectral albedo of
% snow, I, Pure snow. Journal of the Atmospheric Sciences, 37, 2712-2733.
% doi:10.1175/1520-0469(1980)037<2712:AMFTSA>2.0.CO;2
% Warren, S. G., & Wiscombe, W. J. (1980). A model for the spectral albedo of
% snow, II, Snow containing atmospheric aerosols. Journal of the Atmospheric
% Sciences, 37, 2734-2745. doi:10.1175/1520-0469(1980)037<2734:AMFTSA>2.0.CO;2
% Wiscombe, W. J. (1976). The delta-Eddington approximation for a vertically
% inhomogeneous atmosphere (NCAR Technical Note NCAR/TN-121+STR).
% doi:10.5065/D65H7D6Z

% inputs
p = inputParser;
addRequired(p,'mu0',@(x) isempty(x) || (isnumeric(x) && all(x(:)>=0) && all(x(:)<=1)))
addRequired(p,'omega',@(x) isnumeric(x) && all(x(:)>=0) && all(x(:)<=1))
addRequired(p,'g',@(x) isnumeric(x) && all(x(:)>=0) && all(x(:)<=1))
addRequired(p,'method',@ischar)
defaultR0 = 0;
defaultAirmass = [];
defaultTau = Inf;
addParameter(p,'R0',defaultR0,@(x) isnumeric(x) && all(x(:)>=0) && all(x(:)<=1))
addParameter(p,'airmass',defaultAirmass,@(x) isnumeric(x) && all(x(:)>0))
addParameter(p,'tau',defaultTau,@(x) isnumeric(x) && all(x(:)>0))
parse(p,mu0,omega,g,method,varargin{:})

iStruct = parseTwoStreamInput(p.Results);

% choose routine depending on method
if contains(iStruct.method,'hybrid','IgnoreCase',true)
    [Refl,Trans,Beam] = mwHybrid(iStruct);
elseif contains(iStruct.method,'edd','IgnoreCase',true) ||...
        contains(iStruct.method,'ww','IgnoreCase',true) ||...
        contains(iStruct.method,'wisc','IgnoreCase',true) ||...
        contains(iStruct.method,'warr','IgnoreCase',true)
    [Refl,Trans,Beam] = deltaEdd1(iStruct);
else
    error('method ''%s'' not recognized',iStruct.method)
end
if isfield(iStruct,'origSize')
    Refl = reshape(Refl,iStruct.origSize);
    if ~isempty(Trans)
        Trans = reshape(Trans,iStruct.origSize);
    end
    Beam = reshape(Beam,iStruct.origSize);
end
end