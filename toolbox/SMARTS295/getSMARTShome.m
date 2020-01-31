function folder = getSMARTShome()
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
    else
        warning('hostname ''%s'' not recognized, need to add to function %s, but now setting to default',...
            host,mfilename);
        folder = defaultValue;
    end
end

end