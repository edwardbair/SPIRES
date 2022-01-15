function [folder] = getSMARTShome()
% [folder] = getSMARTShome()
%   return name of folder where SMARTS program lives and where input and
%   output files go

vn = getSMARTSversion;

switch vn
    case 295
        folder = getSMARTS295home;
    case 298
        folder = getSMARTS298home;
    otherwise
        if isnumeric(version)
            error('SMARTS version ''%d'' not recognized',version)
        else
            error('SMARTS version ''%s'' not recognized',version)
        end
end
end

function folder = getSMARTS295home()
%returns the home of the SMARTS software
%(depends on host)

defaultValue = 'C:\MyProgramFiles\SMARTS_295_PC';
[status,host] = system('hostname');
if status~=0
    warning('hostname command failed, setting SMARTShome to default')
    folder = defaultValue;
else
    host = strtrim(host);
    if contains(host,'quinaya','IgnoreCase',true)
        folder = 'C:\raid\Program Files\SMARTS_295_PC';
    elseif contains(host,'QUMIAF4','IgnoreCase',true)
        folder = 'C:\MyProgramFiles\SMARTS_295_PC';
    elseif contains(host,'tlun','IgnoreCase',true)
        folder = '/home/snowhydro/nbair/software/SMARTS/SMARTS_295_Linux';
    else
        warning('hostname ''%s'' not recognized, need to add to function %s, but now setting to default',...
            host,mfilename);
        folder = defaultValue;
    end
end
end

function folder = getSMARTS298home()
%returns the home of the SMARTS software
%(depends on host)

defaultValue = 'C:\MyProgramFiles\SMARTS_298_PC';
[status,host] = system('hostname');
if status~=0
    warning('hostname command failed, setting SMARTShome to default')
    folder = defaultValue;
else
    host = strtrim(host);
    if contains(host,'quinaya','IgnoreCase',true)
        folder = 'C:\raid\Program Files\SMARTS_298_PC';
    elseif contains(host,'QUMIAF4','IgnoreCase',true)
        folder = 'C:\MyProgramFiles\SMARTS_298_PC';
    elseif contains(host,'tlun','IgnoreCase',true)
        folder= '/home/snowhydro/nbair/software/SMARTS/SMARTS_295_Linux';
    else
        warning('hostname ''%s'' not recognized, need to add to function %s, but now setting to default',...
            host,mfilename);
        folder = defaultValue;
    end
end
end