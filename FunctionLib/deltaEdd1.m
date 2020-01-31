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