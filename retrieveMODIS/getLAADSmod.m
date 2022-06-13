function [ S ] = getLAADSmod(TOKEN,folder, prefix, date7dig, tile )
%retrieves files from LAADS DAAC website and puts them in specified folder
%
%
% input
%   TOKEN - get from LAADS, may need to be generated w/ EarthData login,
%   i.e.,
%   https://ladsweb.modaps.eosdis.nasa.gov/tools-and-services/data-download-scripts/#generate-token
%   folder - where to put output files
%   prefix - e.g. MOD09GA or MYD09GA
%   date7dig - ISO 7-digit date
%   tile - MODIS tile in form 'hNNvNN'    
%
% output
%   structure with status of commands

%modURL = 'https://e4ftl01.cr.usgs.gov/';
modURL ='https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/';
matdate = iso2datenum(date7dig);
[y,m,d] = datevec(matdate);
doy=datenum([y,m,d])-datenum([y,1,1])+1;

% URL
switch prefix
    case 'MOD09GA'
        url = fullfile(modURL,'61','MOD09GA');
    case 'MOD09GQ'
        url = fullfile(modURL,'61','MOD09GQ');
    otherwise
        error('unrecognized prefix %s',prefix)
end

fullURL = [fullfile(url,num2str(y),num2str(doy,'%03d')),'/'];

str=['wget -v -e robots=off --no-check-certificate --ignore-length --no-use-server-timestamps -r -l Inf -np -nd -R .html,.tmp -nH --cut-dirs=3 '...
    '--header "Authorization: Bearer ' TOKEN '"'...
        ' -P "' folder '" -A "*' ...
        num2str(y,'%04d') num2str(doy,'%03d') '*' tile '*.hdf" ' '"' fullURL '"'];
[s,c] = system(str);
S.statusM = s;
S.cmdoutM = c;

%check integrity
wc=['*' num2str(y,'%04d') num2str(doy,'%03d') '*' tile '*.hdf'];
d=dir(fullfile(folder,wc));

if length(d)==1
    fname=fullfile(d.folder,d.name);
    fprintf('checking integrity of %s\n',fname)
try
    GetMOD09GA(fname,'allbands');
    fprintf('integrity check passed for %s\n',fname)
catch %download again
    warning('integrity check failed, re-downloading %s\n',wc)
    [s,c] = system(str);
    S.statusM = s;
    S.cmdoutM = c;
    try
        GetMOD09GA(fname,'allbands');
    catch
        warning('integrity check failed for 2nd time, deleting %s\n',wc)
        delete(fname);
    end
end

end