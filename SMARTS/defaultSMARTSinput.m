function prescription = defaultSMARTSinput(RefAtmos,varargin)
% prescription = defaultSMARTSinput(RefAtmos,varargin)
%defaultSMARTSinput - generate input for SMARTSMain with default values
%
%Required input
%   RefAtmos - reference atmosphere118 to use, options are (case-insensitive):
%       'mlw' - mid-latitude winter
%       'sas' - Sub-Arctic summer
%       others to be added
%Optional input, name/value pairs
%   'cosZ' - cosine of solar zenith on horizontal surface
%   'altit' - elevation of surface, km (default 3000 for 'mlw', 0 for 'sas', 0 for 'ussa')
%   'height' - elevation above surface, km (default 0)
%
%For other options, use the more general routine SetSMARTS

p = inputParser;
defaultZen = 48.19;
defaultCosZ = cosd(defaultZen);
if contains(RefAtmos,'mlw','IgnoreCase',true)
    defaultALTIT = 3000;
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
addParameter(p,'wlmn',280,@isnumeric)
addParameter(p,'wlmx',4000,@isnumeric)
addParameter(p,'albedo',[],@(x) isscalar(x) && x>=0 && x<=1)
parse(p,RefAtmos,varargin{:});

switch lower(p.Results.RefAtmos)
    case 'mlw'
        if isempty(p.Results.albedo)
            argc = {'COMNT','mid-latitude_winter','ISPR',2,'LATIT',37,...
                'ALTIT',p.Results.altit,'LATIT',37.6,'IATMOS',1,'IH2O',1,'IALBDX',30,...
                'IPRT',2,'IMASS',0,'ZENIT',acosd(p.Results.cosz),'ATMOS','MLW',...
                'HEIGHT',p.Results.height,'WLMN',p.Results.wlmn,'WLMX',p.Results.wlmx};
        else
            argc = {'COMNT','mid-latitude_winter','ISPR',2,'LATIT',37,...
                'ALTIT',p.Results.altit,'LATIT',37.6,'IATMOS',1,'IH2O',1,'IALBDX',-1,...
                'RHOX',p.Results.albedo,...
                'IPRT',2,'IMASS',0,'ZENIT',acosd(p.Results.cosz),'ATMOS','MLW',...
                'HEIGHT',p.Results.height,...
                'WLMN',p.Results.wlmn,'WLMX',p.Results.wlmx};
        end
    case 'sas'
        z0 = 0;
        T0 = 273;
        P0 = 1013.25;
        spr = AirPressure(P0,T0,z0,[],p.Results.altit); %#ok<NASGU>
        argc = {'COMNT','sub-Arctic_summer','ISPR',2,'LATIT',75,...
            'LATIT',70,'IATMOS',1,'IH2O',1,'IALBDX',0,...
            'IPRT',2,'IMASS',0,'ZENIT',acosd(p.Results.cosz),'ATMOS','SAS',...
            'ALTIT',p.Results.altit,'HEIGHT',p.Results.height,...
            'WLMN',p.Results.wlmn,'WLMX',p.Results.wlmx};
    case 'ussa'
        z0 = 0;
        T0 = 288;
        P0 = 1013.25;
        spr = AirPressure(P0,T0,z0,[],p.Results.altit);
        argc = {'COMNT','US_Standard_Atmosphere','ISPR',1,'SPR',spr,...
            'LATIT',38,'IATMOS',1,'IH2O',1,'ISPCTR',0,'IALBDX',0,...
            'IPRT',2,'IMASS',0,'ZENIT',acosd(p.Results.cosz),'ATMOS','USSA',...
            'ALTIT',p.Results.altit,'HEIGHT',p.Results.height,...
            'WLMN',p.Results.wlmn,'WLMX',p.Results.wlmx};
    otherwise
        error('RefAtmos ''%s'' unknown',p.Results.RefAtmos)
end

%set the prescription
prescription = SetSMARTS(argc{:});
end