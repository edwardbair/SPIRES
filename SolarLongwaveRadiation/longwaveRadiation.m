function [ L ] = longwaveRadiation( TK, solarIndex, varargin )
% [ L ] = longwaveRadiation( TK, solarIndex, {'precipWater', w} or {'vaporPressure',e} )
%compute incoming longwave radiation based on methods of Dilley-O'Brien, Prata,
% and Crawford-Duchon
%
% Input
%   TK - air temperature, deg K
%   solarIndex - (actual / clear sky)
% And, either name-value pair (case-insensitive)
%   'vaporPressure' in hPa
%   or 'precipWater' in cm
%
% Output
%   longwave radiation - W/m^2
%
% Calculate clear sky longwave, if w is available, using the formula in
% Dilley & O'Brien, 1998, QJRMS doi: 10.1002/qj.49712454903.
%
% If instead vapor pressure is available, then Prata, 1996, QJRMS
% doi: 10.1002/qj.49712253306.
%
% Then adjust based on the solar index s Crawford & Duchon, 1999, JAM
% doi: 10.1175/1520-0450(1999)038<0474:AIPFEE>2.0.CO;2.
%
% This strategy is based on the recommendation in Flerchinger et al., 2009,
% WRR doi: 10.1029/2008WR007394.

p = inputParser;
V = {lower(varargin{1}) varargin{2}};
addRequired(p,'TK',@(x) isnumeric(x) && all((x(:)>100)))
addRequired(p,'solarIndex',@isnumeric);
addParameter(p,'precipwater',[],@isnumeric);
addParameter(p,'vaporpressure',[],@isnumeric);
parse(p,TK,solarIndex,V{:});

assert(xor(isempty(p.Results.precipwater),isempty(p.Results.vaporpressure)),...
    'either precipitable water or vapor pressure must be specified, but not both')

S = radiationConstants;
sigma = S.StefanBoltzmann;
Tmelt = 273.15;

% w = precipitable water in cm
% s = solar index

% clear values, if w is available - eq 16 and Table 1 in Dilley-O'Brien
if ~isempty(p.Results.precipwater)
    xi = p.Results.precipwater;
else
    % or if vapor pressure is available - eq 20 in Prata
    vaporPressure = p.Results.vaporpressure;
    xi = 46.5*(vaporPressure./TK);
end
Lclr = 59.38+113.7.*(TK/Tmelt).^6+96.96*sqrt(xi/25);
emissClear = Lclr./(sigma*TK.^4);

% adjust for solar index - eqn 3 to 11 in Crawford-Duchon
emiss = (1-solarIndex)+solarIndex.*emissClear;

L = sigma*emiss.*TK.^4;

end