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
    if any(isnan(Refl))
        t = isnan(Refl);
        iStruct.semiInfinite = true;
        [rx,tx,bx] = mwHybrid(iStruct);
        Refl(t) = rx(t);
        Trans(t) = tx(t);
        Beam(t) = bx(t);
    end
elseif contains(iStruct.method,'edd','IgnoreCase',true) ||...
        contains(iStruct.method,'ww','IgnoreCase',true) ||...
        contains(iStruct.method,'wisc','IgnoreCase',true) ||...
        contains(iStruct.method,'warr','IgnoreCase',true)
    [Refl,Trans,Beam] = deltaEdd1(iStruct);
    if any(isnan(Refl))
        t = isnan(Refl);
        iStruct.semiInfinite = true;
        [rx,tx,bx] = deltaEdd1(iStruct);
        Refl(t) = rx(t);
        Trans(t) = tx(t);
        Beam(t) = bx(t);
    end
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

function [refl,trans,beam] = deltaEdd1(S)
% [refl,trans,beam] = deltaEdd1(S)
%Wiscombe-Warren (1980) for 1-layer delta Eddington
%Input
% S - structure from parseTwostreamInputs
%Output
% reflectance, diffuse transmittance, beam transmittance

%intermediate and transformed variables
gstar = S.g./(1+S.g);
ostar = (1-S.g.^2).*S.omega./(1-S.g.^2.*S.omega);
astar = 1-ostar.*gstar;
bstar = gstar./astar;
xisq = 3*astar.*(1-ostar);
xi = sqrt(xisq);
P = 2*xi./(3.*astar);

if S.semiInfinite
    if S.Direct % WW80, eq 4
        refl = (ostar./(1+P)).*((1-bstar.*xi.*S.mu0)./(1+xi.*S.mu0));
    else % WW80, eq 7
        refl = (2*ostar./(1+P)).*(((1+bstar)./xisq).*(xi-log(1+xi))-bstar/2);
    end
    trans = zeros(size(refl));
    beam = zeros(size(refl));
else
    % more intermediate variables
    tstar = (1-S.omega.*S.g).*S.tau0;
    gam = (1-S.R0)./(1+S.R0);
    Qp = (gam+P).*exp(xi.*tstar);
    Qm = (gam-P).*exp(-xi.*tstar);
    Q = (1+P).*Qp-(1-P).*Qm;
    if S.Direct % WW80, eq 3
        if isfield(S,'pathLength')
            pathLength = S.pathLength;
        else
            pathLength = 1./S.mu0;
        end
        nRef = 2*(P.*(1-gam+ostar.*bstar)+ostar.*(1+bstar).*...
            ((gam.*xi.*S.mu0-P)./(1-xi.^2.*S.mu0.^2))).*...
            exp(-tstar.*pathLength)-ostar.*bstar.*(Qp-Qm)+...
            ostar.*(1+bstar).*(Qp./(1+xi.*S.mu0)-Qm./(1-xi.*S.mu0));
        beam = exp(-tstar.*pathLength);
        % hack-hack, haven't figured this out yet so using M-W hybrid model
        % with the delta-Eddington transformation of omega, g, tau (seems
        % to work for the reflectance too but we stay with delta-Eddington)
        [~,trans,~] = twostream(S.mu0,ostar,gstar,'hybrid','tau',tstar,'R0',S.R0);
        refl = nRef./Q;
    else % WW80, integrating eq 3 rather than using WW80 ex 6
        %I was not able to get Wiscombe-Warren eq 6 to work, but instead I
        %integrated (2*mu0*albedo(mu0) dmu0 from 0 to 1) in a way
        %that is only valid for xi<1 (difficult to integrate eq 3
        %because some xi>1 and (1-mu0^2*xi^2) is in the denominator
        %For values where xi>=1, I simply use Wiscombe's observation that
        %diffuse properties are similar to direct at cosines between about
        %0.6 and 0.7
        
        t = xi<1;
        nRef = (ostar(t).*xi(t).*(2*(1+bstar(t)).*(Qm(t)+Qp(t))+bstar(t).*(Qm(t)-Qp(t)).*xi(t))+...
            2*exp(-tstar(t)).*xi(t).*(-2*(1+bstar(t)).*gam(t).*ostar(t)+(-1+gam(t)-bstar(t).*ostar(t)).*P(t).*(-1+tstar(t)).*xi(t))+...
            2*(2*(1+bstar(t)).*ostar(t).*P(t)+2*(1+bstar(t)).*gam(t).*ostar(t).*tstar(t).*xi(t)+...
            (1-gam(t)+bstar(t).*ostar(t)).*P(t).*tstar(t).^2.*xi(t).^2).*expint(tstar(t))+...
            2*(1+bstar(t)).*exp(-tstar(t).*xi(t)).*ostar(t).*(gam(t)-P(t)).*expint(tstar(t).*(1-xi(t)))+...
            2*(1+bstar(t)).*ostar(t).*(-(exp(tstar(t).*xi(t)).*(gam(t)+P(t)).*expint(tstar(t).*(1+xi(t))))+Qm(t).*log(1-xi(t))-...
            Qp(t).*log(1+xi(t))))./xi(t).^2;
        
        % hack here, diffuse properties very similar to direct at cos=0.65;
        [refl,trans,~] = twostream(0.65,ostar,gstar,'hybrid','tau',tstar,'R0',S.R0);
        beam = zeros(size(refl));
        % replace some refl values (for xi<1) with the analytic delta-Eddington integral
        refl(t) = nRef./Q(t);
    end
end
end

function [refl,trans,beam] = mwHybrid(S)
% [refl,trans,beam] = mwHybrid(S)
%Meador-Weaver hybrid method for calculating reflectance, diffuse
%transmittance, and direct transmittance of a single layer

%Input
% S - structure from parseTwostreamInputs
%Output
% reflectance, diffuse transmittance, and direct transmittance

assert(S.Direct,'diffuse reflectance not available from mwHybrid (yet) because difficult to analytically integrate')

gam = mwgamma(S.mu0,S.omega,S.g);
xi = sqrt((gam(:,1)-gam(:,2)).*(gam(:,1)+gam(:,2)));
if S.semiInfinite
    refl = S.omega.*(gam(:,3).*(xi+gam(:,1)-gam(:,2))+gam(:,2))./...
        ((xi+gam(:,1)).*(1+xi.*S.mu0));
    trans = zeros(size(refl));
    beam = zeros(size(refl));
    % corrections for conservative scattering
    t = S.omega==1;
    if any(t)
        refl(t) = 1;
    end
    t = S.omega==0;
    if any(t)
        refl(t) = 0;
    end
else
    % more intermediate variables
    if isfield(S,'pathLength')
        pathLength = S.pathLength;
    else
        pathLength = 1./S.mu0;
    end
    alph1 = gam(:,1).*gam(:,4) + gam(:,2).*gam(:,3);
    alph2 = gam(:,2).*gam(:,4) + gam(:,1).*gam(:,3);
    em = exp(-S.tau0.*xi);
    et = exp(-S.tau0.*pathLength);
    ep = exp(S.tau0.*xi);
    gpx = xi + gam(:,1);
    opx = S.mu0 .* xi + 1;
    omx = 1 - S.mu0 .* xi;
    gmx = gam(:,1) - xi;
    rm = gam(:,2) - gmx.*S.R0;
    rp = gam(:,2) - gpx.*S.R0;
    
    % denominator for reflectance and transmittance
    denrt = ep.*gpx.*rm - em.*gmx.*rp;
    
    % reflectance
    refl = (S.omega.*(ep.*rm.*(gam(:,3).*xi+alph2)./opx -...
        em.*rp.*(alph2-gam(:,3).*xi)./omx) + 2*et.*gam(:,2).*...
        (S.R0-((alph1.*S.R0-alph2).*S.mu0+gam(:,4).*S.R0+gam(:,3)).*S.omega./...
        (omx.*opx)).*xi)./denrt;
    
    % transmittance
    trans = (et.*(ep.*gpx.*(gam(:,2)-S.omega.*(alph2-gam(:,3).*xi)./omx)-...
        em.*gmx.*(gam(:,2)-S.omega.*(gam(:,3).*xi+alph2)./opx))+...
        2*gam(:,2).*(alph1.*S.mu0+gam(:,4)).*S.omega.*xi./(omx.*opx))./denrt;
    trans(trans>1) = 1-8*eps;
    
    % direct transmittance
    beam = et;
end
end

function [ gam ] = mwgamma( mu0, omega, g)
% [ gam ] = mwgamma( mu0, omega, g )
%mwgamma - Meador-Weaver gamma values
%
% Input
%   mu0 - cosine of illumination angle
%   omega - single-scattering albedo
%   g - asymmetry factor
%
% Output
%   gamma - values as a column vectors (i.e. Nx4 matrix)

% check input ranges
assert(isequal(size(mu0),size(omega),size(g)),...
    'input mu0, omega, g must be same size')
assert(all(mu0>=0 & mu0<=1), 'mu0 must be 0-1')
assert(all(omega>=0 & omega<=1), 'omega must be 0-1')
assert(all(g>=0 & g<=1), 'g=%g, must be 0-1')

gam = zeros(length(g),4);

% Meador-Weaver Table 1, some in Horner polynomial form
b0 = betanaught(mu0, g);
hd = 4+g.^2.*(-4+4*mu0); % denominator for both gam(1:2), = 4*(1-g^2*(1-mu0))
gam(:,1) = (7-4*omega+g.*(-3*omega+g.*(-3+4*b0.*omega+3*g.*omega)))./hd;
gam(:,2) = (-1+4*omega+g.*(-3*omega+g.*(1+(-4+4*b0).*omega+3*g.*omega)))./hd;
gam(:,3) = b0;
gam(:,4) = 1-gam(:,3);

% insure gam(:1)=gam(:2) for conservative scattering (i.e. no rounding
% error)
t = omega==1;
if any(t)
    gam(t,2) = gam(t,1);
end
end

function beta0 = betanaught( mu0, g )
% beta0 = betanaught( mu0, g )
%compute Meador-Weaver beta_0 value
%   Equations 3 and 4
%
% Input
%   mu0 - cosine of illumination angle
%   g - scattering asymmetry parameter

% parameters
MAXNO = 2048;
TOL = eps('single')*10;

% check inputs
assert(all(mu0 >= 0 & mu0 <= 1))
assert(all(g >= 0 & g <= 1))

% sum until convergence; we use the even terms only for the recursive
% calculation of the Legendre polynomials
bsum = zeros(size(g));
for k=1:length(g)
    % Legendre polynomials of degree 0 and 1
    pnm2 = 1;
    pnm1 = mu0(k);
    % first coefficients and initial sum
    fm = -1/8;
    gn = 7 * g(k)^3;
    bsum(k) = 3 * g(k) * mu0(k) / 2;
    if g(k) ~= 0 && mu0(k) ~= 0
        for n=2:MAXNO
            % order n Legendre polynomial
            pn = ((2 * n - 1) * mu0(k) * pnm1 + (1 - n) * pnm2) / n;
            if mod(n,2) == 1 % odd n
                last = bsum(k);
                bsum(k) = last + gn * fm * pn;
                if abs((bsum(k) - last) / bsum(k)) < TOL && bsum(k) <= 1
                    break;
                end
                % recursively find next f(m) and gn coefficients n = 2 * m + 1
                m = (n - 1) / 2;
                fm = fm * (-(2 * m + 1) / (2 * (m + 2)));
                gn = gn * (g(k)^2 * (4 * m + 7) / (4 * m + 3));
            end
            % ready to compute next Legendre polynomial
            pnm2 = pnm1;
            pnm1 = pn;
        end
        if n >= MAXNO
            conv = (bsum(k)-last)/bsum(k);
            if conv>TOL*10
                warning('%s: %s - mu0=%g g=%g bsum=%g last=%g conv=%g fm=%g',...
                    'betanaught', 'no convergence',...
                    mu0(k), g(k), bsum(k), last, conv, fm);
            end
            if (bsum(k) > 1)
                bsum(k) = 1;
            end
        else
            assert(bsum(k) >= 0 && bsum(k) <= 1,'numerical error');
        end
    else
        bsum(k) = 0;
    end
end

beta0 = (1-bsum)/2;

end

function S = parseTwoStreamInput(R)
% S = parseTwoStreamInput(R)
% parse inputs for the twostream codes, values returned in structure S

% angles, omega, and g must be same size
% direct or diffuse reflectance?
if isempty(R.mu0)
    S.Direct = false;
    [S.omega,S.g] = checkSizes(R.omega,R.g);
else
    S.Direct = true;
    [S.mu0,S.omega,S.g] = checkSizes(R.mu0,R.omega,R.g);
end

% semi-infinite?
if (isscalar(R.tau) && isinf(R.tau)) ||...
        (~isscalar(R.tau) && all(isinf(R.tau(:))))
    S.semiInfinite = true;
else
    S.semiInfinite = false;
    [S.tau0,S.R0,~] = checkSizes(R.tau,R.R0,S.omega);
end
S.method = R.method;
if ~isempty(R.airmass)
    assert(~isempty(R.mu0),'if ''airmass'' is specified, so must cosine of incidence angle')
    [S.pathLength,~] = checkSizes(R.airmass,S.omega);
end

% if inputs are matrices, convert to vectors but save original size
if size(S.omega,1)~=1 && size(S.omega,2) ~=1
    S.origSize = size(S.omega);
    fn = fieldnames(S);
    for k=1:length(fn)
        if ~(islogical(S.(fn{k})) || ischar(S.(fn{k})) ||...
                contains(fn{k},'size','IgnoreCase',true))
            S.(fn{k}) = S.(fn{k})(:);
        end
    end
end
end
    