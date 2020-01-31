function prescription = defaultSMARTSinput(RefAtmos,varargin)
% prescription = defaultSMARTSinput(RefAtmos,varargin)
%defaultSMARTSinput - generate input for SMARTS295Main with default values
%
%Required input
%   RefAtmos - reference atmosphere to use, options are (case-insensitive):
%       'mlw' - mid-latitude winter
%       'sas' - Sub-Arctic summer
%       others to be added
%Optional input, name/value pairs
%   'cosZ' - cosine of solar zenith on horizontal surface
%   'cosS' - cosine of slope for tilted surface
%   'altit' - elevation of surface, km (default 3 for 'mlw', 0 for 'sas', 0 for 'ussa')
%   'height' - elevation above surface, km (default 0)
%
%For other options, use the more general routine SetSMARTS295

p = inputParser;
defaultZen = 48.19;
defaultCosZ = cosd(defaultZen);
if contains(RefAtmos,'mlw','IgnoreCase',true)
    defaultALTIT = 3;
else
    defaultALTIT = 0;
end
defaultHEIGHT = 0;
addRequired(p,'RefAtmos',@ischar)
addParameter(p,'cosz',defaultCosZ,...
    @(x) isscalar(x) && isnumeric(x) && x>=0 && x<=1)
addParameter(p,validatestring('altit',{'alt','elev','altit'}),...
    defaultALTIT,@(x) isscalar(x) && isnumeric(x) && x<=70)
addParameter(p,validatestring('height',{'ht','height'}),...
    defaultHEIGHT,@(x) isscalar(x) && isnumeric(x) && x<=70)
parse(p,RefAtmos,varargin{:});

ioutValues = [1 2 3 5 6 7 21 30];

switch lower(p.Results.RefAtmos)
    case 'mlw'
        P0 = 700;
        z0 = 3000;
        T0 = 265;
        spr = AirPressure(P0,T0,z0,[],p.Results.altit*1000);
        argc = {'COMNT','mid-latitude winter','ISPR',1,'SPR',spr,...
            'ALTIT',p.Results.altit,'LATIT',37.6,'IATMOS',1,'IH2O',1,'ISPCTR',0,'IALBDX',30,...
            'IPRT',2,'IMASS',0,'ZENIT',acosd(p.Results.cosz),'ATMOS','MLW','IOUT',...
            ioutValues,'HEIGHT',p.Results.height};
    case 'sas'
        z0 = 0;
        T0 = 273;
        P0 = 1013.25;
        spr = AirPressure(P0,T0,z0,[],p.Results.altit*1000);
        argc = {'COMNT','sub-Arctic summer','ISPR',1,'SPR',spr,...
            'LATIT',70,'IATMOS',1,'IH2O',1,'ISPCTR',0,'IALBDX',0,...
            'IPRT',2,'IMASS',0,'ZENIT',acosd(p.Results.cosz),'ATMOS','SAS','IOUT',...
            ioutValues,'ALTIT',p.Results.altit,'HEIGHT',p.Results.height};
    case 'ussa'
        z0 = 0;
        T0 = 288;
        P0 = 1013.25;
        spr = AirPressure(P0,T0,z0,[],p.Results.altit*1000);
        argc = {'COMNT','US Standard Atmosphere','ISPR',1,'SPR',spr,...
            'LATIT',38,'IATMOS',1,'IH2O',1,'ISPCTR',0,'IALBDX',0,...
            'IPRT',2,'IMASS',0,'ZENIT',acosd(p.Results.cosz),'ATMOS','USSA','IOUT',...
            ioutValues,'ALTIT',p.Results.altit,'HEIGHT',p.Results.height};
    otherwise
        error('RefAtmos ''%s'' unknown',p.Results.RefAtmos)
end

%set the prescription
prescription = SetSMARTS295(argc{:});
end