function S = radiationConstants()
% S = radiationConstants
%returns the radiation constants widely used in the physics of
%electromagnetic radiation:
% Planck
% SpeedOfLight
% Boltzmann
% StefanBoltzmann

narginchk(0,0)
nargoutchk(0,1)

% fundamental constants:
h = 6.62607015e-34; % Planck, J s
k = 1.380649e-23; % Boltzmann, J/K
c = 299792458; % speed of light, m/s

% Stefan-Boltzmann, based on integral of Planck equation from 0 to Inf
sigma = (2*k^4*pi^5)/(15*c^2*h^3);

S.Description = 'fundamental radiation constants';
S.Planck = h;
S.Boltzmann = k;
S.SpeedOfLight = c;
S.StefanBoltzmann = sigma;

end