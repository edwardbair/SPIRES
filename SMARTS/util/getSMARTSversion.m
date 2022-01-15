function SMARTS_version = getSMARTSversion()
% SMARTS_version = getSMARTSversion()
%   get SMARTS version from global variable
global SESSION_SMARTS_version

defaultVersion = 298;
if isempty(SESSION_SMARTS_version)
    warning('SMARTS version not set, run SetSMARTSversion, defaulting to ''%d''',...
        defaultVersion)
    SMARTS_version = defaultVersion;
else
    SMARTS_version = SESSION_SMARTS_version;
end
end
