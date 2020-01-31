function [ H, L ] = BusingerDyer( zu, zt, z0, T, T0, e, e0, u, u0, P )
% [ H, L ] = BusingerDyer( zu, zt, z0, T, T0, e, e0, u, u0, P )
%sensible and latent heat exchange from the Businger-Dyer formulation
%
% Input - scalars or vectors/matrices, but if not scalars must be same size
%   zu height of upper wind speed measurement (m)
%   zt height of temp & vapor pressure measurements (m)
%   z0 either height of lower wind measurement, or roughness length (m)
%   T temperature at zt (K)
%   T0 temperature at z0 (K), surface temperature if z0 is roughness length
%   e vapor pressure at zt (kPa)
%   e0 vapor pressure at z0 (kPa), can be at surface
%   u wind speed at zu (m/s)
%   u0 wind speed at lower height, or zero if z0 is roughness length
%   P air pressure (kPa)
%
% Output
%   H sensible heat exchange, W/m^2
%   L latent heat exchange, W/m^2
%
% Reference: Fleagle & Businger, An Introduction to Atmospheric Physics,
% 2nd ed., 1980

% anonymous functions
Tmelt = 273.16;
Rstar = 8.31432; % J/mole/deg
mol_H2O = 0.01801528; % kg/mole
mol_Air = 0.0289644; % kg/mole
lh_vap = @(T) (2.5e6 - 2.95573e3 *(T - Tmelt)); % J/kg
spec_hum = @(e,P) ((e).*mol_H2O./(mol_Air.*(P)+(e)*(mol_H2O-mol_Air)));
vir_temp = @(T,e,P) (T./(1-(1-mol_H2O/mol_Air).*(e./P)));
gas_den = @(P,m,T) (P.*m./(Rstar*T));

% check for errors in input
assert(all(zu(:)>z0(:)), 'zu must be greater than z0')
assert(all(zt(:)>z0(:)), 'zt must be greater than z0')
assert(all(T(:)>200), 'T must be in Kelvin')
assert(all(T0(:)>200), 'T0 must be in Kelvin')

% sizes of input variables - insert check here

% correct vapor pressures above saturation
satvpT = SaturationVaporPressure(T,'water');
tooWet = e>satvpT;
if any(tooWet(:))
    e(tooWet) = satvpT(tooWet);
    warning('some vapor pressures reset to saturation vapor pressure')
end
ice = T0<=Tmelt;
satvp0(~ice) = SaturationVaporPressure(T0(~ice),'water');
satvp0(ice) = SaturationVaporPressure(T0(ice),'ice');
tooWet = e0>satvp0;
if any(tooWet(:))
    e0(tooWet) = satvp0(tooWet);
    warning('some surface vapor pressures reset to saturation vapor pressure')
end

% geometric mean T, e
tgm = sqrt(T.*T0);
egm = sqrt(e.*e0);

% specific heat of air
cp_air = 1005; % J kg^-1 deg^-1
cp = cp_air;

% latent heat of vaporization, temperature dependent
xlh = lh_vap(tgm);

% specific humidity, dimensionless
q = spec_hum(e,P);
q0 = spec_hum(e0,P);

% air density
aden = gas_den(100*P,mol_Air,vir_temp(tgm,egm,P));

[ustar,Tstar,qstar] = bdinit();

H = cp*aden.*ustar.*Tstar;
L = aden.*xlh.*qstar.*ustar;

% nested function for Businger-Dyer values
    function [ustar,Tstar,qstar] = bdinit()
        k = 0.41; % von Karman constant
        gravity = 9.80665;
        dalr = gravity/cp_air;
        alpha = 16;
        beta = 4.7;
        RiCrit = 1 / beta;
        lzz0 = log(zu./z0);
        
        % convert T to potential temperature, wrt surface
        potT = T+dalr*zt;
        
        % bulk Richardson number at geometric mean height
        Ri = RichardsonNumber(zt,z0,potT,T0,u,u0);
        
        % Businger-Dyer psi functions
        psi1 = zeros(size(Ri));
        psi2 = zeros(size(Ri));
        
        % 4 flow regimes
        laminar = Ri<=RiCrit;
        zeroRi = Ri==0;
        negRi = Ri<0;
        posRi = ~laminar & Ri>0;
        
        % negative Richardson number
        zeta = Ri(negRi);
        % x = 1-phi_sub_m
        y = sqrt(1-alpha*zeta);
        x = sqrt(y);
        % the Businger-Dyer psi-functions
        psi1(negRi) = 2.*log((1+x)/2)+log((1+y)/2)-2*atan(x)+pi/2;
        psi2(negRi) = 2*log((1+y)/2);
        
        % positive Richardson number
        zeta = Ri(posRi)./(1-beta*Ri(posRi));
        psi1(posRi) = -beta*zeta;
        psi2(posRi) = psi1(posRi);
        
        % output values
        Tstar = k*(T-T0)./(lzz0-psi2);
        ustar = k*(u-u0)./(lzz0-psi1);
        qstar = k*(q-q0)./(lzz0-psi2);
        % all zero if laminar
        Tstar(laminar|zeroRi) = 0;
        ustar(laminar) = 0;
        qstar(laminar) = 0;
        
    end
end

function Ri = RichardsonNumber(z2,z1,t2,t1,u2,u1)
% input
%   z2, z1 - heights of temperature and wind
%   t2, t1 - virtual temperatures at z2, z1
%   u2, u1 - winds at z2, z1 (u1 can be zero)

assert(all(z2(:)>z1(:)) & all(z2(:)>0),'some upper heights > lower heights')
assert(all(u2(:)>u1(:)) & all(u2(:)>0) & all(u1(:)>=0),...
    'some upper winds > lower')

% geometric mean virtual temperature
tgm = sqrt(t2.*t1);

gravity = 9.80665;
Ri = gravity*(z2-z1).*(t2-t1)./(tgm.*(u2-u1).^2);
end