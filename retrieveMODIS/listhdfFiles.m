function hdffiles = listhdfFiles(folder,prefix,thisYear,tile)
fhdf = dir(fullfile(folder,[prefix '*A' num2str(thisYear) '*' tile '*hdf']));
if isempty(fhdf)
    warning('no files in %s matching ''%s''',folder,...
        [prefix '*A' num2str(thisYear) '*' tile '*hdf'])
    warning('will download files as part of processing')
    hdffiles = '';
else
    hdffiles = cell(length(fhdf),1);
    for k=1:length(hdffiles)
        hdffiles{k} = fullfile(folder,fhdf(k).name);
    end
end
end