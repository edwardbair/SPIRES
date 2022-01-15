  function C=SetSMARTS(varargin)
% C= SetSMARTS298(Name/Value pairs)
%Creates structure for input file and then uses struct2SMARTSinput to print
%the input file smarts298.inp.txt, which can be run with smarts298.exe
%Input options are described in detail in the SMARTS manualiout
%smarts298-users-manual-pc.pdf, in the Documentation folder of the SMARTS_298_PC root folder.
%The following brief descriptions use the same variable names as the documentation.
%
%Name/Value pairs, in any order, most CARD# are required with exceptions as noted, CARD#a are optional
%CARD1
% 'COMNT' - brief comment about this run, <= 64 characters including spaces
%CARD2
% 'ISPR' - 0, 1, or 2 to indicate how surface pressure is specified (see docs)
%CARD2a
% 'SPR' - surface pressure in mb, default 1013.25, ISPR = 0 or 1
% 'ALTIT' - surface altitude in m, default 0, ISPR = 1 or 2 (note that the
%   value passed to SMARTS is in km)
% 'HEIGHT' - height above surface altitude, m, default 0, ISPR = 1 or 2
%   (note that the value passed to SMARTS is in km)
% 'LATIT' - latitude, degrees, default 45, ISPR = 2
%CARD3 atmosphere
% 'IATMOS' - 0 for user-selected atmosphere, 1 to use a reference atmosphere
%CARD3a if 'IATMOS' = 0
% 'TAIR' - air temperature (C) at reference height
% 'RH' - relative humidity (%)
% 'SEASON' - either 'WINTER' or 'SUMMER', default 'WINTER', case-insensitive
% 'TDAY' - mean daily temperature (C)
%CARD3a if 'IATMOS' = 1
% 'ATMOS' - name of selected atmosphere, 4 char max, default 'MLW' (see documentation for descriptions)
%CARD4 precipitable water
% 'IH2O' - 0 if precipitable water specified, 1 if from 'ATMOS'
% (if precipitable water is available, should specify 'IATMOS' as 0, the CARD3a must have TAIR etc)
%CARD4a if 'IH2O' = 0
% 'W' - precipitable water in cm, max 12
%CARD5 ozone, use default as options not implemented yet
%CARD6 other gases, use default as options not implemented yet
%CARD7 CO2
% 'qCO2' - columnar concentation, ppmv, can use default 400
%CARD7a extraterrestrial spectrum
% 'ISPCTR' - default 9 for version 2.9.8, default 0 for 2.9.5,
%   see documentation for options -1 to 9, here briefly
%   0 - Gueymard 2004
%   1 - Gueymard unpublished
%   2-7 - various MODTRAN options
%   8 - ASTM E490, 2000
%   9 - Gueymard 2018ab
%CARD8 aerosol model
% 'AEROS' - see documentation for options, default is 'S&F_RURAL' (Shuttle & Fine Rural)
%CARD9 turbidity
% 'ITURB' - see documentation for values from 0 to 5, default 0
%CARD9a turbidity values see documentation for possibilities, default for ITURB=0 is TAU5=0.084
%values corresponding to 'ITURB' are:
%0 'TAU5', 1 'BETA', 2 'BCHUEP', 3 'RANGE', 4 'VISI', 5 'TAU550'
%CARD10 regional (10 km) albedo
% 'IALBDX' - values -1 to 66, -1 a user-defined gray albedo specified by
%   'RHOX', 0 or +1 a user-defined spectral albedo in file Albedo.dat,
%   2 to 66 see albedo tables in documentation
%CARD10a - value of 'RHOX' 0 to 1
%CARD10b - 'ITILT', 0 if not tilted surface, 1 if tilted
%CARD10c - if 'ITILT' = 1, then values for 'IALBDG', 'TILT' 0 to 90 or -999
%   for sun-tracking photometer, and 'WAZIM' orientation 0 to 360 clockwise
%CARD10d - sames as CARD10a except for 'RHOG'
%CARD11 = values of 'WLMN' minimum wavelength in nm default 280,
%   'WLMX' maximum wavelength default 4000,
%   'SUNCOR' Earth-Sun radius vector default 1,
%   'SOLARC' solar constant default 1361.1 W/m^2
%CARD12 - output options for variable 'IPRT'
%   0 - only integrated results are printed to File smarts298.out.txt
%   1 - integrated and spectral results are printed to File smarts298.out.txt
%   2 - spectral results are printed to separate File smarts298.ext.txt
%CARD12a - values of WPMN, WPMX, INTVL
%   WPMN = minimum wavelength, must be >= WLMN
%   WPMX - maximum wavelength, must be <= WLMX
%   INTVL - wavelength interval, can be coarser than SMARTS native values
%       but if finer set to actual intervals, as is case when default 0.5 is
%       used (usually best to leave this alone)
%CARD12b - auto-calculated
%CARD12c - 'IOUT' variables to output, vector with numbers from 1 to 43, see
%   documentation for names that correspond to values
%CARD13 - currently defaulted, for circumsolar calculations
%CARD14 - currently defaulted, for specific spectral instruments
%CARD15 - currently defaulted, for photosynthetically active radiation
%CARD16 - currently defaulted, for UV radiation
%CARD17 - 'IMASS' airmass option
%   0, 'ZENIT' and 'AZIM' on CARD17a
%   1, 'ELEV' and 'AZIM' on CARD17a
%   2, 'AMASS' on CARD17a
%   3, 'YEAR', 'MONTH', 'DAY', 'HOUR', 'LATIT', 'LONGIT', 'ZONE' on CARD17a
%   4, 'MONTH', 'LATIT', 'DSTEP' for daily calculation on CARD17a

version = getSMARTSversion;
switch version
    case 295
        defaultSolarSpectrum = 0;
        upperLimit = 8;
    case 298
        defaultSolarSpectrum = 9;
        upperLimit = 9;
end
       

p = inputParser;
% CARD1
addParameter(p,'COMNT','default_comment',@(x) ischar(x) && length(x)<=64)
% CARD2
addParameter(p,'ISPR',[],@(x) isscalar(x) && x>=0 && x<=2 && x==round(x))
% CARD2a, depending on value in CARD2
defaultSPR = 1013.25;
addParameter(p,'SPR',defaultSPR,@(x) isscalar(x) && x>0)  % ISPR 0 or 1
addParameter(p,'ALTIT',0,@(x) isscalar(x) && x<=100000) % ISPR 1 or 2
addParameter(p,'HEIGHT',0,@(x) isscalar(x) && x<=100000) % ISPR 1 or 2
addParameter(p,'LATIT',45,@(x) isscalar(x) && x>=-90 && x<=90)
%CARD3
addParameter(p,'IATMOS',[],@(x) isscalar(x) && (x==1 || x==0))
%CARD3a
addParameter(p,'TAIR',[],@(x) isempty(x) || (isscalar(x) && x>=-120 && x<=50))
addParameter(p,'RH',[],@(x) isempty(x) || (isscalar(x) && x>=0 && x<=100))
addParameter(p,'SEASON','WINTER',@(x) ischar(x) && (strcmpi(x,'WINTER') ||...
    strcmpi(x,'SUMMER')))
addParameter(p,'TDAY',[],@(x) isempty(x) || (isscalar(x) && x>=-120 && x<=50))
addParameter(p,'ATMOS','USSA',@(x) ischar(x) && length(x)<=4)
%CARD4
addParameter(p,'IH2O',1,@(x) isscalar(x) && x==0 || x==1) % doc recommends against IH2O=1, so not allowed
%CARD4a
addParameter(p,'W',[],@(x) isempty(x) || (isscalar(x) && x>=0 && x<=12))
%CARD5 [options for ozone not yet implemented]
addParameter(p,'IO3',1,@(x) isscalar(x) && isnumeric(x) && x==0 || x==1)
%CARD6 [options for other gases not yet implemented]
addParameter(p,'IGAS',1,@(x) isscalar(x) && isnumeric(x) && x==0 || x==1)
%CARD7
addParameter(p,'qCO2',400,@(x) isscalar(x) && isnumeric(x))
%CARD7a
addParameter(p,'ISPCTR',defaultSolarSpectrum,@(x) isscalar(x) && isnumeric(x)...
    && x>=-1 && x<=upperLimit)
%CARD8 aerosol model
addParameter(p,'AEROS','S&F_RURAL',@ischar)
%CARD9 turbidity
addParameter(p,'ITURB',0,@(x) isnumeric(x) && x>=0 && x<=5)
addParameter(p,'TAU5',0.084,@(x) isnumeric(x) && x>0)
addParameter(p,'BETA',[],@isnumeric)
addParameter(p,'BSCHUEP',[],@isnumeric)
addParameter(p,'RANGE',[],@isnumeric)
addParameter(p,'VISI',[],@isnumeric)
addParameter(p,'TAU550',[],@isnumeric)
%CARD10 albedo
addParameter(p,'IALBDX',[],@(x) isempty(x) || (isnumeric(x) && x>=-1 && x<=66))
addParameter(p,'RHOX',.5,@(x) isnumeric(x) && x>=0 && x<=1)
addParameter(p,'RHOG',.5,@(x) isnumeric(x) && x>=0 && x<=1)
addParameter(p,'ITILT',0,@(x) isnumeric(x) && (x==0 || x==1))
addParameter(p,'IALBDG',[],@(x) isempty(x) || (isnumeric(x) && x>=-1 && x<=66))
addParameter(p,'TILT',[],@(x) isempty(x) ||...
    (isnumeric(x) && (x==-999 || (x>=0 && x<=90))))
addParameter(p,'WAZIM',[],@(x) isempty(x) || (isnumeric(x) && x>=0 && x<=360))
%CARD11 wavelengths and solar information
addParameter(p,'WLMN',280,@(x) isnumeric(x) && x>=280)
addParameter(p,'WLMX',4000,@(x) isnumeric(x) && x<=4000)
addParameter(p,'SUNCOR',1,@(x) isnumeric(x) && x>=0.96 && x<=1.04)
addParameter(p,'SOLARC',1361.1,@(x) isnumeric(x) && x>0)
%CARD12 output options
addParameter(p,'IPRT',0,@(x) isnumeric(x) && x>=0 && x<=3 && x==round(x))
%CARD12a
addParameter(p,'WPMN',[],@isnumeric);
addParameter(p,'WPMX',[],@isnumeric);
addParameter(p,'INTVL',0.5,@isnumeric);
%CARD12c
addParameter(p,'IOUT',[1:8 21 30 32 34],@(x) isnumeric(x) && all(x>=1) && all(x<=43))
%CARD12b - IOTOT is calculated so this field is normally empty, just
%provided for completeness
addParameter(p,'IOTOT',[],@isnumeric);
%CARD13 through CARD16 defaulted
addParameter(p,'ICIRC',0,@isnumeric)
addParameter(p,'ISCAN',0,@isnumeric)
addParameter(p,'ILLUM',0,@isnumeric)
addParameter(p,'IUV',0,@isnumeric)
%CARD17
addParameter(p,'IMASS',[],@(x) isnumeric(x) && x>=0 && x<=4 && x==round(x))
addParameter(p,'ZENIT',45,@(x) isnumeric(x) && x>=0 && x<=90)
addParameter(p,'AZIM',180,@(x) isnumeric(x) && x>=0 && x<=360)
addParameter(p,'ELEV',[],@(x) isempty(x) || (isnumeric(x) && x>=0 && x<=90))
addParameter(p,'AMASS',[],@(x) isempty(x) || (isnumeric(x) && x>=0))
% note: 'LATIT' covered possibly for CARD2a
addParameter(p,'YEAR',[],@(x) isempty(x) || (isnumeric(x) && x==round(x)))
addParameter(p,'MONTH',[],@(x) isempty(x) ||...
    (isnumeric(x) && x==round(x) && x>=1 && x<=12))
addParameter(p,'DAY',[],@(x) isempty(x) ||...
    (isnumeric(x) && x==round(x) && x>=1 && x<=31))
addParameter(p,'HOUR',[],@(x) isempty(x) ||...
    (isnumeric(x) && x>=0 && x<=24))
addParameter(p,'LONGIT',0,@(x) isscalar(x) && x>=-180 && x<=180)
addParameter(p,'ZONE',0,@(x) isscalar(x) && x>=-12 && x<=12)
addParameter(p,'DSTEP',15,@(x) isscalar(x) && x>=5 && x<=60 && mod(60,x)==0)

parse(p,varargin{:})

% Note: docs say replace blanks w '_' but apparent not necessary
CARD1.Value = p.Results.COMNT; %#ok<*STRNU>
CARD1.Name = 'COMNT';

CARD2.Value = int8(p.Results.ISPR);
CARD2.Name = 'ISPR';
assert(~isempty(CARD2.Value),'''ISPR'' must be specified, either 0, 1, or 2')
switch CARD2.Value
    case 0
        card2a.Value = p.Results.SPR;
        card2a.Name = 'SPR';
    case 1
        card2a.Value = [p.Results.SPR p.Results.ALTIT/1000 p.Results.HEIGHT/1000];
        % check consistency of pressure and elevation
        if p.Results.SPR~=defaultSPR && p.Results.ALTIT==0
            warning('check values of ''SPR'' and ''ALTIT'', possibly inconsistent')
        end
        assert(sum(card2a.Value(2:3))<=100,'''ALTIT''+''HEIGHT'' must be <=100 km')
        card2a.Name = 'SPR ALTIT HEIGHT';
    case 2
        card2a.Value = [p.Results.LATIT p.Results.ALTIT/1000 p.Results.HEIGHT/1000];
        assert(sum(card2a.Value(2:3))<=100,'''ALTIT''+''HEIGHT'' must be <=100 km')
        card2a.Name = 'LATIT ALTIT HEIGHT';
end

CARD3.Value = int8(p.Results.IATMOS);
CARD3.Name = 'IATMOS';
assert(~isempty(CARD3.Value),'''IATMOS'' must be specified, either 0 or 1')
switch CARD3.Value
    case 0
        % cell vector instead of numeric vector because ischar(SEASON)
        card3a.Value = {p.Results.TAIR,p.Results.RH,p.Results.SEASON,p.Results.TDAY};
        card3a.Name = 'TAIR RH SEASON TDAY';
    case 1
        card3a.Value = validatestring(p.Results.ATMOS,...
            {'AS','AW','MLS','MLW','SAS','SAW','STS','STW','TRL','USSA'});
        card3a.Name = 'ATMOS';
end

CARD4.Value = int8(p.Results.IH2O);
CARD4.Name = 'IH2O';
switch CARD4.Value
    case 0
        card4a.Value = p.Results.W;
        assert(~isempty(p.Results.W),'if ''IH2O''=0, then ''W'' must be specified')
        card4a.Name = 'W';
    case 1
        if ~isempty(p.Results.W)
            warning('Because ''IH2O''=1, input value of ''W'' is ignored, value from ''ATMOS'' is used')
        end
end

CARD5.Value = int8(p.Results.IO3);
CARD5.Name = 'IO3';
assert(CARD5.Value==1,...
    '''%s'' values other than 1 not yet implemented, use 0 for ''%s'' and read documentation',...
    CARD5.Name,CARD5.Name)

CARD6.Value = int8(p.Results.IGAS);
CARD6.Name = 'IGAS';
assert(CARD6.Value==1,...
    '''%s'' values other than 1 not yet implemented, use 0 for ''%s'' and read documentation',...
    CARD6.Name,CARD6.Name)

CARD7.Value = p.Results.qCO2;
CARD7.Name = 'qCO2';
card7a.Value = int8(p.Results.ISPCTR);
card7a.Name = 'ISPCTR';

CARD8.Name = 'AEROS';
CARD8.Value = validatestring(p.Results.AEROS,...
    {'B&D_C','B&D_C1','DESERT_MAX','DESERT_MIN','S&F_MARIT','S&F_RURAL',...
    'S&F_TROPO','S&F_URBAN','SRA_CONTL','SRA_MARIT','SRA_URBAN','USER'});
assert(~strcmpi(CARD8.Value,'USER'),...
    'option for ''USER'' defined aerosol not yet implemented')

CARD9.Value = int8(p.Results.ITURB);
CARD9.Name = 'ITURB';
switch CARD9.Value
    case 0
        card9a.Value = p.Results.TAU5;
        card9a.Name = 'TAU5';
    case 1
        if isempty(p.Results.BETA)
            card9a.Value = p.Results.TAU5;
            CARD9.Value = 0;
            card9a.Name = 'TAU5';
        else
            card9a.Value = p.Results.BETA;
            card9a.Name = 'BETA';
        end
    case 2
        if isempty(p.Results.BSCHUEP)
            card9a.Value = p.Results.TAU5;
            CARD9.Value = 0;
            card9a.Name = 'TAU5';
        else
            card9a.Value = p.Results.BSCHUEP;
            card9a.Name = 'BSCHUEP';
        end
    case 3
        if isempty(p.Results.RANGE)
            card9a.Value = p.Results.TAU5;
            CARD9.Value = 0;
            card9a.Name = 'TAU5';
        else
            card9a.Value = p.Results.RANGE;
            card9a.Name = 'RANGE';
        end
    case 4
        if isempty(p.Results.VISI)
            CARD9.Value = 0;
            card9a.Value = p.Results.TAU5;
            card9a.Name = 'TAU5';
        else
            card9a.Value = p.Results.VISI;
            card9a.Name = 'VISI';
        end
    case 5
        if iempty(p.Results.TAU550)
            card9a.Value = p.Results.TAU5;
            CARD9.Value = 0;
            card9a.Name = 'TAU5';
        else
            card9a.Value = p.Results.TAU550;
            card9a.Name = 'RAU5';
        end
end

CARD10.Value = int8(p.Results.IALBDX);
assert(~isempty(p.Results.IALBDX),'must enter value for ''IALBDX''')
CARD10.Name = 'IALBDX';
if CARD10.Value==-1
    card10a.Value = p.Results.RHOX;
    card10a.Name = 'RHOX';
end
card10b.Value = int8(p.Results.ITILT);
card10b.Name = 'ITILT';
if card10b.Value==1
    card10c.Value = [p.Results.IALBDG p.Results.TILT p.Results.WAZIM];
    card10c.Name = 'ALBDG TILT WAZIM';
    card10d.Value = p.Results.RHOG;
    card10d.Name = 'RHOG';
end

assert(p.Results.WLMN<p.Results.WLMX,...
    'min wavelength must be < max wavelength -- check ''WLMN'' and ''WLMX''')
CARD11.Value = [p.Results.WLMN p.Results.WLMX p.Results.SUNCOR p.Results.SOLARC];
CARD11.Name = 'WLMN WLMX SUNCOR SOLARC';

CARD12.Value = int8(p.Results.IPRT);
CARD12.Name = 'IPRT';
% using defaults for CARD12a
if CARD12.Value>=1
    if isempty(p.Results.WPMN)
        wpmn = p.Results.WLMN;
    else
        wpmn = max([p.Results.WLMN p.Results.WPMN]);
    end
    if isempty(p.Results.WPMX)
        wpmx = p.Results.WLMX;
    else
        wpmx = min([p.Results.WLMX p.Results.WPMX]);
    end
    card12a.Value = [wpmn wpmx p.Results.INTVL];
    card12a.Name = 'WLMN WLMX INTVL';
end
if CARD12.Value>=2
    card12c.Value = p.Results.IOUT;
    card12c.Name = 'IOUT';
    card12b.Value = length(card12c.Value);
    card12b.Name = 'IOTOT';
end

%CARD13 through CARD16 currently defaulted
CARD13.Value = int8(p.Results.ICIRC);
CARD14.Value = int8(p.Results.ISCAN);
CARD15.Value = int8(p.Results.ILLUM);
CARD16.Value = int8(p.Results.IUV);
CARD13.Name = 'ICIRC';
CARD14.Name = 'ISCAN';
CARD15.Name = 'ILLUM';
CARD16.Name = 'IUV';

CARD17.Value = int8(p.Results.IMASS);
CARD17.Name = 'IMASS';
switch CARD17.Value
    case 0
        card17a.Value = [p.Results.ZENIT p.Results.AZIM];
        card17a.Name = 'ZENIT AZIM';
    case 1
        card17a.Value = [p.Results.ELEV p.Results.AZIM];
        card17a.Name = 'ELEV AZIM';
    case 2
        card17a.Value = p.Results.AMASS;
        card17a.Name = 'AMASS';
    case 3
        card17a.Value = [p.Results.YEAR p.Results.MONTH p.Results.DAY...
            p.Results.HOUR p.Results.LATIT p.Results.LONGIT p.Results.ZONE];
        card17a.Name = 'YEAR MONTH DAY HOUR LATIT LONGIT ZONE';
    case 4
        card17a.Value = [p.Results.MONTH p.Results.LATIT. p.Results.DSTEP];
        card17a.Name = 'MONTH LATIT DSTEP';
end

% now populate the structure with all the "cards"
N = 17;
L = {'a','b','c','d'};
for n=1:N
    cname = ['CARD' num2str(n)];
    if exist(cname,'var')==1
        C.(cname) = eval(cname);
        for k=1:length(L)
            suppname = ['card' num2str(n) L{k}];
            if exist(suppname,'var')==1
                C.(suppname) = eval(suppname);
            end
        end
    end
end
end