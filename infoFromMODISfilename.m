function fileTable = infoFromMODISfilename( f )
% fileTable = infoFromMODISfilename( f )
%
%   Extract datenum (and maybe time), prefix, and tile (if gridded) from
%   MODIS-like filenames
%
% Input
%   f - cell vector of filenames, can be full path
%
% Output - table variable with following columns, sorted by matdate
%   matdate - list of dates in MATLAB datenum format
%   possibly enddate, if 2 dates in file name
%   file - list of files
%   prefix = MODIS product suite
%   tile - MODIS tile names, if gridded

assert(iscell(f)||ischar(f),...
    'input filelist must be cell vector or character string');
if ischar(f)
    f = {f};
end
N = length(f);
dm = zeros(N,length(dateFromName(f{1}))); % possibly a date range?
T = cell(N,1);
P = cell(N,1);
if isrow(f) && N>1
    f = f';
end

haveTile = false;
for k=1:N
    dm(k,:) = dateFromName(f{k});
    tn = tileFromName(f{k});
    if ~isempty(char(tn))
        haveTile = true;
    end
    T(k) = tn;
    P(k) = prefixFromName(f{k});
end

% create table, sort in chronological order, and return
if haveTile % tiles available
    if size(dm,2)==1
        Tbl = table(dm,P,T,f,'VariableNames',{'matdate','prefix','tile','file'});
    else
        Tbl = table(dm(:,1),dm(:,2),P,T,f,'VariableNames',...
            {'matdate','enddate','prefix','tile','file'});
    end
else % no tiles
    if size(dm,2)==1
        Tbl = table(dm,P,f,'VariableNames',{'matdate','prefix','file'});
    else
        Tbl = table(dm(:,1),dm(:,2),P,f,'VariableNames',...
            {'matdate','enddate','prefix','file'});
    end
end
    fileTable = sortrows(Tbl);
end

% removes folder information, through final slash
    function s=fnstart(c)
        k = strfind(c,'\');
        if isempty(k)
            k = strfind(c,'/');
        end
        if isempty(k)
            s = c;
        else
            s = c(k(end)+1:end);
        end
    end

% date from file name, which is in form .AYYYYDDD or form _AYYYYDDD or form
% _AYYYYDDD_YYYYDDD, or if time of day AYYYYDDD.HHMM
    function D=dateFromName(c)
        c = fnstart(c);
        k = strfindi(c,'.A');
        if isempty(k)
            k = strfindi(c,'_A');
        end
        if isempty(k)
            error('date not found in file name %s (in .Ayyyyddd or _Ayyyyddd format)',c)
        end
        kk = k(1)+2;
        yyyyddd = str2double(c(kk:kk+6));
        D = iso2datenum(yyyyddd);
        % does a second date follow? (i.e. does this file cover a date range)?
        nexts = c(kk+7:end);
        if strcmp(nexts(1),'_')
            % are the next 7 characters digits
            possibleDate = nexts(2:8);
            if all(isstrprop(possibleDate,'digit'))
                yyyyddd = str2double(possibleDate);
                D(2) = iso2datenum(yyyyddd);
            end
        else % or does a time of day follow?
            if strcmp(nexts(1),'.')
                possibletime = nexts(2:5);
                if all(isstrprop(possibletime,'digit'))
                    hhmm = str2double(possibletime);
                    D = D+floor(hhmm/100)/24+mod(hhmm,100)/(24*60);
                end
            end
        end
    end

% tile from MODIS file name
    function T=tileFromName(c)
        c = fnstart(c);
        k = strfindi(c,'_h');
        if isempty(k)
            k = strfindi(c,'.h');
            if isempty(k) || strcmp(c(k(1):end),'.hdf') % no tile in this file name
                T = {''};
                return
            end
        end
        kk = k(1)+1;
        T = {c(kk:kk+5)};
    end

% prefix from MODIS file name
    function P=prefixFromName(c)
        c = fnstart(c);
        k = strfind(c,'.');
        if isempty(k) || k(1)>15
            k = strfind(c,'_');
        end
        P = {c(1:k-1)};
    end