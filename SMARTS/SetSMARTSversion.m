function SetSMARTSversion(version)
% SetSMARTSversion(version)
%   set SMARTS version for all runs this MATLAB session unless cleared

global SESSION_SMARTS_version

p = inputParser;
addRequired(p,'version',@(x) ischar(x) || isscalar(x) &&...
    (isstring(x) || isnumeric(x)))
parse(p,version)
version = p.Results.version;

if ischar(version) || isstring(version)
    vn = str2double(replace(version,'.',''));
else
    vn = version;
end

SESSION_SMARTS_version = uint16(vn);

end

