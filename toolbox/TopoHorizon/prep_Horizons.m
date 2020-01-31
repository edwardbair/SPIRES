function prep_Horizons(demfile, nNodes, varargin )
% prep_Horizons(demfile, nNodes )
%   Consolidate the results from HorizonsAtDirections, for additional
%   parallel processing
%
% Input
%   demfile - full path name to DEM file
%   nNodes - number of nodes used (to decide how to parse for parallel
%       processing)
%
% Output
%   written to files in same folder as demfile
%   Each contains original horizon header and piece of horizon grid
%

if nargin<2
    if isdeployed
        disp(['usage: ' mfilename ' demfile nNodes'])
    else
        disp(['usage: ' mfilename '(demfile, nNodes)'])
    end
    return
end

p = inputParser;
addRequired(p,'demfile',@ischar)
addRequired(p,'nNodes',@(x) ischar(x) || isnumeric(x))
parse(p,demfile,nNodes)

demfile = p.Results.demfile;
[dempath,demname,ext]= fileparts(demfile);
if isempty(dempath)
    dempath = '.';
end

% number of nodes to prepare for
if isnumeric(p.Results.nNodes)
    N = p.Results.nNodes;
else
    N = str2double(p.Results.nNodes);
end
nNodes = N;

% input horizon folder and output diary file
folder = dempath;
if isdeployed
diaryFolder = fullfile(folder,'output');
    if exist(diaryFolder,'file')~=7
        mkdir(diaryFolder);
    end
    StartAzureDiary(mfilename,diaryFolder,demname,...
        ['nNodes' num2str(nNodes,'%03d')]);
end

% input files
f = dir(fullfile(folder,[demname '*H0' '*.mat']));
for k=1:length(f)
    horzfile{k,1} = f(k).name; %#ok<AGROW>
end
nFiles = length(horzfile);
disp('horizon files')
disp(horzfile)
disp(['elevation file ' fullfile(dempath,[demname ext])])

% slopes and aspects
ZS = load(demfile);
S = SlopeAzimuth(ZS);

% load the first file and initialize the output matrices
Horz = load(fullfile(folder,horzfile{1}));
sH = zeros(size(Horz.sinH,1),size(Horz.sinH,2),nFiles,'like',Horz.sinH);
Az = zeros(size(Horz.azm,1),size(Horz.azm,2),nFiles,'like',Horz.azm);
sH(:,:,1) = Horz.sinH;
Az(:,:,1) = Horz.azm;

% rest of the files
for k=2:nFiles
    H = load(fullfile(folder,horzfile{k}));
    sH(:,:,k) = H.sinH;
    Az(:,:,k) = H.azm;
end
assert(isequal(size(sH),size(Az)),'sH and Az cubes must be same size')
if nnz(sH>1)
    warning('%d horizon sines are greater than 1.0',nnz(sH>1))
end

origN = size(sH);
% change sH & Az to size [nHorz nRows*nCols], i.e. reshape then transpose
sH = reshape(sH,origN(1)*origN(2),origN(3))';
Az = reshape(Az,origN(1)*origN(2),origN(3))';

% reshape slope & aspect into linear vectors
XS = reshape(S.Slope,1,numel(S.Slope));
XA = reshape(S.Aspect,1,numel(S.Aspect));

% note - no need to sort, as postprocess_Horizons does that

% break up Az & sH into separate files for parallel processing
N = size(sH);
nodeUsage = 2; % number of times each node gets used
siz = bestblk(N,round(N(2)/(nNodes*nodeUsage)));
nOut = ceil(N(2)/siz(2));
s = zeros(nOut,1);
s(1) = 1;
for k=2:nOut
    s(k) = s(k-1)+siz(2);
end
e = s+siz(2)-1;
e(e>N(2)) = N(2);

% write the files
fname = [demname '_Hblock'];
for k=1:nOut
    Horz.sinH = sH(:,s(k):e(k));
    Horz.azm = Az(:,s(k):e(k));
    Horz.slope = XS(:,s(k):e(k));
    Horz.aspect = XA(:,s(k):e(k));
    SaveBlock(folder,[fname num2str(k,'%04d')],Horz)
end

end

function SaveBlock(folder,filename,Horz )
% SaveBlock(folder, filename,Horz)
% writes temporary horizon block information
% modify this function depending on where data are to be stored
%
% folder - name of folder
% filename - name of file (no folder info)
% Horz - horizon structure

save(fullfile(folder,filename),'-struct','Horz','-v7.3')

end