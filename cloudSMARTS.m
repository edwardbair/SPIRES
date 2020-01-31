function outStruct = cloudSMARTS(smartsHome,outputFolder,varargin)
% outStruct = cloudSMARTS(smartsHome,outputFolder,name/value pairs)
%
%cloudSMARTS runs SMARTS295 with a cloud layer
%
%Inputs to this function, same as for SMARTS295Main but with a cloud layer
%added, with ALTIT set to the ground elevation and HEIGHT the height of the
%cloud above that, both in km
%The first two arguments:
% smartsHome - the folder where SMARTS is installed
% outputFolder - where results are saved as .mat files, can be empty
%   (either '' or []) to just get the results in the output structure
%   without writing them to a file.
%Third and following arguments
% Optionally, the 3rd argument can be an input prescription structure for
%   SMARTS, perhaps modfied from the output of a previous run of SMARTS.
%   If not, the 3rd and following arguments consist of the input cell
%   vector to run SMARTS295Main, then followed by a cell vector whos first
%   entry is 'cloudType', followed by the type of cloud ('watercloud',
%   'icecloud', or 'mixedCloud') and then the cloud properties
% Optionally, if the 3rd argument is a prescription structure for SMARTS,
% then the 4th argument can be a prescription structure for
% SnowCloudSpectralRefl, or the 4th and following arguments can start with
% 'cloudType', followed by cloud properties

%% parse the variable argument cell vector
useSMARTSvarargin = true;
useCloudvarargin = true;
if isstruct(varargin{1})
    smartsP = varargin{1};
    useSMARTSvarargin = false;
    varargin = varargin(2:end);
    if isstruct(varargin{1})
        cloudP = varargin{1};
        useCloudvarargin = false;
    end
end
if useCloudvarargin
    % separate the cloud varargin cell vector
    t = strcmpi(varargin,'cloudType');
    assert(any(t),'variable argument vector must include ''cloudType''')
    k = find(t);
    assert(length(k)==1,'variable argument vector must have just one instance of ''cloudType''')
    varCloud = varargin(k+1:end);
    if useSMARTSvarargin
        varSMARTS = varargin(1:k-1);
    end
end

% set the SMARTS prescriptions if not in the calling arguments
if useSMARTSvarargin
    smartsP = SetSMARTS295(varSMARTS{:});
end

%% run SMARTS once to get the current prescription
sS = SMARTS295Main(smartsHome,outputFolder,smartsP);
argc = argScriptFromStructure(sS.prescription);
% set the cloud prescription, with cosZ from the arguments
if useCloudvarargin
    holdv = cat(2,varCloud,{'wavelength',sS.spectralTbl.waveL,...
        'waveUnit','nm'});
    cloudP1 = SetSnowCloud(holdv{:});
    if ~isfield(cloudP1,'cosZ')
        t = strcmpi(argc,'zenit');
        if any(t)
            k = find(t);
            solZ = argc{k+1};
        end
        cloudP1.cosZ = cosd(solZ);
    end
    holdv = cat(2,varCloud,{'wavelength',sS.spectralTbl.waveL,...
        'waveUnit','nm','cosZ',[]});
    cloudP2 = SetSnowCloud(holdv{:});
else
    cloudP.wavelength = sS.spectralTbl.waveL;
    cloudP.waveUnit = 'nm';
    cloudP1 = cloudP;
    cloudP2 = cloudP;
    if ~isfield(cloudP1,'cosZ')
        t = strcmpi(argc,'zenit');
        if any(t)
            k = find(t);
            solZ = argc{k+1};
        end
        cloudP1.cosZ = cosd(solZ);
    end
    cloudP2.cosZ = [];
end

% now set variables needed for input, consistent with SnowCloudSpectralRefl
% and SMARTS solar zenith
t = strcmpi(argc,'zenit');
if any(t)
    k = find(t);
    argc{k+1} = acosd(cloudP1.cosZ);
else
    argc = cat(2,argc,{'ZENIT',acosd(cloudP1.cosZ)});
end
% input variables must include HEIGHT
t = strcmpi(argc,'height');
if any(t)
    k = find(t);
    if argc{k+1}==0
        argc{k+1} = 1; % 1 km default cloud height above surface
    end
else
    argc = cat(2,argc,{'HEIGHT',1});
end
% output variables, must include HorzDiffuse, HorzDirect, and ReflZone
t = strcmpi(argc,'iout');
if any(t)
    k = find(t);
    argc{k+1} = unique(cat(2,argc{k+1},[3 5 30]));
    t = strcmpi(argc,'iotot');
    m = find(t);
    argc{m+1} = numel(argc{k+1});
else
    argc = cat(2,argc,{'IOUT',[3 5 30],'IOTOT',3});
end

% and run SMARTS
sS = SMARTS295Main(smartsHome,outputFolder,argc{:});


%% run SnowCloudSpectralRefl in both direct and diffuse mode
[Sdir,Pdir] = SnowCloudSpectralRefl(cloudP1);
Sdif = SnowCloudSpectralRefl(cloudP2);

% irradiance under the cloud
Rdir = sS.spectralTbl.HorzDirect.*Sdir.beam;
Rdif = sS.spectralTbl.HorzDirect.*Sdir.trans +...
    sS.spectralTbl.HorzDiffuse.*Sdif.trans;

sS.spectralTblUnderCloud = table(sS.spectralTbl.waveL,Rdir,Rdif,...
    'VariableNames',{'waveL','HorzDirect','HorzDiffuse'});
sS.cloudPrescription = Pdir;
outStruct = sS;
end
