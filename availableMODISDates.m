function [date7dig,files] = availableMODISDates(fileList,firstday,lastday)
fT = infoFromMODISfilename(fileList);
matdate = fT{:,'matdate'};
files = fT{:,'file'};
date7dig = datenum2iso(unique(matdate),7);
dayvec = rem(date7dig,1000);
if ~isempty(lastday)
    t = dayvec>lastday;
    date7dig(t) = [];
end
if ~isempty(firstday)
    t = dayvec<firstday;
    date7dig(t) = [];
end
if isempty(date7dig)
    files = '';
end
end