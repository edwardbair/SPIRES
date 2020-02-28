function [ fileTable] = infoFromLandsatFilename( f )
% [ fileTable] = infoFromLandsatFilename( f )
%
%   Extract datenum, filetype, sensor id, and WRS-2 path/row from Landsat-like filenames
%
% Input
%   f - cell vector of filenames, can be full path
%
%
% Output table colums, rows sorted by matdate
%   matdate - sorted list of dates in MATLAB datenum format
%   filetype - Can be LS8 bands 1-11, quality band, MTL metadata,
%       or cloudmask shapefile(as long as shapfile is named appropriately).
%       In same order as matdate, or a folder, or .gz or .tar.gz
%   WRS2 - cell vector of WRS-2 path/row scene identifiers in same order
%       as matdate
%   sortedFiles - list of files in same order as matdate
%
%   extension - filetype extension (.tif, .txt, .shp, ext..)
%   satID - cell vector of Landsat sensor & satellite in same order as matdate
%   sceneID - scene ID for surrent scene.
%
%
% Example sceneID - LC80290302015343LGN00
% Example productID - LC08_L1GT_029030_20151209_20160131_01_RT
%
%
% january 2018 updates:
% can handle filenames without extensions
% can handle classic pre-collection landsat sceneID's and the new landsat
% collection ID scehma.

assert(iscell(f)||ischar(f),...
    'input filelist must be cell vector or character string');
if ischar(f)
    f = {f};
end

%empty vectors for outputs
N = length(f);
sceneIDidx = false(N,1);
productIDidx = false(N,1);
matdate = zeros(N,1);
fileType = cell(N,1);
extension = cell(N,1);
wrs2 = cell(N,1);
satID = cell(N,1);
typeID = cell(N,1);

if isrow(f) && N>1
    f = f';
end


%determine if Scene ID or Landsat Product Identifer file naming convention



for k=1:N
    
    
    %Check that file name is valid
    c = f{k};
    c = fnstart(c);
    
    %1st 10 alphabetic and numeric elements in correct location
    sceneIDcorrect = [1 1 0 0 0 0 0 0 0 0];
    productIDcorrect = [1 1 0 0 0 1 0 1 1 0];
    
    tf = isstrprop(c(1:10),'alpha');
    
    sceneIDidx(k) = all(tf == sceneIDcorrect);
    productIDidx(k) = all(tf == productIDcorrect);
    assert(sceneIDidx(k) | productIDidx(k), 'incorrect filename format %s', c)
    
    %1st element is L
    assert(c(1) == 'L', 'incorrect landsat file name %s', c)
    
    
    if productIDidx(k)
        assert(str2double(c(4)) == 8, 'file not from landsat 8, from landsat %s', c(3));
        sceneID = {c(1:40)};
        pIDswitch = true;
        typeID(k) = {'productID'};
        
    elseif sceneIDidx(k)
        assert(str2double(c(3)) == 8, 'file not from landsat 8, from landsat %s', c(3));
        sceneID = {c(1:21)};
        pIDswitch = false;
        typeID(k) = {'sceneID'};
    end
    
    
    %output
    
    matdate(k) = dateFromName(f{k},pIDswitch);
    [fileType(k), extension(k)] = filetypeFromName(f{k});
    wrs2(k) = WRS2fromName(f{k},pIDswitch);
    satID(k) = idFromName(f{k},pIDswitch);

end




% create table, sort in chronological order, and return
if length(sceneID)<length(satID)
    sceneID = repmat(sceneID,size(satID));
end
tbl = table(matdate, fileType, wrs2, f, extension, satID, typeID, sceneID, ...
    'VariableNames',{'matdate' 'filetype' 'path_row'...
    'file' 'file_extension' 'satID','typeID', 'sceneID'});
%fileTable = sortrows(tbl);
fileTable = tbl;
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

% date from file name, which is in form YYYYDDD
function D=dateFromName(c,pIDswitch)
c = fnstart(c);

if ~pIDswitch
    yyyyddd = str2double(c(10:16));
    D = iso2datenum(yyyyddd);
else
    yyyymmdd = str2double(c(18:25));
    D = iso2datenum(yyyymmdd);
end


end

% WRS2 path/row from file name, which is in form PPPRRR
function pathrow = WRS2fromName(c,pIDswitch)

if ~pIDswitch
    c = fnstart(c);
    pathrow = {['p' c(4:6) 'r' c(7:9)]};
else
    c = fnstart(c);
    pathrow = {['p' c(11:13) 'r' c(14:16)]};
end


end

% id from Landsat file name
function W=idFromName(c,pIDswitch)

if ~pIDswitch
    c = fnstart(c);
    W = {c(1:3)};
else
    c = fnstart(c);
    W = {c(1:4)};
end


end

% filetype from Landsat file name -NEED TO BE REWRITTEN TO WORK WITH
% PRODUCT ID.
function [ft, extension] =filetypeFromName(f)
c = fnstart(f);
k = strfind(c,'_');
j = strfind(c,'.');

%folder or file? - if folder, must exist.
if ~isdir(f)
    if ~isempty(j)
        extension = {c(j(end)+1:end)};
        if strcmp(extension,'tar')
            ft = {'tar-wad'};
        elseif strcmp(extension,'gz')
            ft = {'compressed tar-wad'};
        elseif ~isempty(k)
            ft = {c(k(end)+1:j(end)-1)};
        end
%     elseif ~isempty(k)
%         ft = {c(k+1:end)};
%         extension = {'unknown'};
    else
        ft = {'nameOnly'};
        extension = {'none'};
        %warning('unreadable filename %s', c)
    end
else
    ft = {'folder'};
    extension = {'none'};
end
end