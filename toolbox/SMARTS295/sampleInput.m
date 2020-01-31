%Examples of ways to craft input for smarts295Main.
%You can adapt these code segments to create an input structure with
%SetSMARTS295.m and then pass the structure to smarts295Main (or
%SnowCloudIntgRefl) or you can pass the cell vector directly to
%smarts295Main)

%Example for CUES site, mid-latitude winter atmosphere, 3 km elevation,
%solar zenith 48.19 (atmos path=1.5), rural aerosols which is the default
% NOTE: THE NAME/VALUE PAIRS CAN BE IN ANY ORDER
argCUES = {'COMNT','test run for SMARTS @ CUES','ISPR',1,'SPR',700,...
    'ALTIT',3,'LATIT',37.6,'IATMOS',1,'IH2O',1,'ISPCTR',0,'IALBDX',30,...
    'IPRT',2,'IMASS',0,'ZENIT',48.19,'ATMOS','MLW','IOUT',...
    [1 2 3 5 21 30]};
% then you can either
Scues1 = SMARTS295Main(getSMARTShome,[],argCUES{:});
% or
Pcues = SetSMARTS295(argCUES{:});
Scues2 = SMARTS295Main(getSMARTShome,[],Pcues);
% results should be the same
disp('Are Scues1 and Scues2 the same?')
isequal(Scues1,Scues2)

%Example for Sub-Arctic summer atmosphere, sea level (the default)
%solar zenith 48.19
argArctic = {'COMNT','test run for SMARTS in Sub-Arctic','ISPR',1,...
    'LATIT',70,'IATMOS',1,'IH2O',1,'ISPCTR',0,'IALBDX',0,...
    'IPRT',2,'IMASS',0,'ZENIT',48.19,'ATMOS','SAS'};
PA = SetSMARTS295(argArctic{:});
SA = SMARTS295Main(getSMARTShome,[],PA);

%you can also recover the whole input argument string without using any
%default values
argNoDefaults = argScriptFromStructure(PA);

%test of cloudSMARTS
argCloud = {'cloudType','watercloud','radius',4,'watereq',1};
