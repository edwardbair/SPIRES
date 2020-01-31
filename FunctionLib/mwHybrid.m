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