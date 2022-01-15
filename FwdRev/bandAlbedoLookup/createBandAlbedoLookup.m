function [F,S,filename] = createBandAlbedoLookup(inputStruct,folder,type,varargin)
% [F,S,file] = createBandAlbedoLookup(inputStruct,folder,type [, name-value pairs])
%   generate lookup functions for each band using output from
%   prepBandAlbedoLookupClean or prepBandAlbedoLookupDirty
%   and append cell vectors for the functions into an output file in the
%   current folder
%
%Input
%   inputStruct - structure with independent variables, band passes, and band
%       IDs
%   folder - where 'checkClean' or 'checkDirt' files were written by
%       prepBandAlbedoLookup{Clean,Dirty}
%   type - 'clean' or 'dirty'
%Optional input, name-value pairs
%   'method' - interpolation method, default 'linear'
%   'extrapolation' - default 'linear' on clean, 'nearest' for dirty
%   'sensor' - default ''
%Output
%   F - cell vector of functions, one for each band
%   S - output structure with band information and about sensor if exist

if strcmpi(type,'clean')
    defaultExtrap = 'linear';
else
    defaultExtrap = 'nearest';
end
p = inputParser;
addRequired(p,'inputStruct',@isstruct)
addRequired(p,'folder',@ischar)
addRequired(p,'type',@(x) ischar(x) && (strcmpi(x,'clean') || strcmpi(x,'dirty')))
addParameter(p,'method','linear',@ischar)
addParameter(p,validatestring('extrapolation',{'extr','extrapolation'}),...
    defaultExtrap,@ischar)
addParameter(p,'sensor','',@ischar)
parse(p,inputStruct,folder,type,varargin{:});

%parse inputStruct
if isfield(inputStruct,'bandPass')
    S.bandPass = inputStruct.bandPass;
    inputStruct = rmfield(inputStruct,'bandPass');
end
if isfield(inputStruct,'bands')
    S.bands = inputStruct.bands;
    inputStruct = rmfield(inputStruct,'bands');
end
if ~isempty(p.Results.sensor)
    S.sensor = p.Results.sensor;
end

switch lower(p.Results.type)
    case 'clean'
        fn = {'elevation','mu0','SSA','muS'};
    case 'dirty'
        fn = {'elevation','mu0','SSA','muS','dust','soot'};
end
nF = length(fn);
S.variables = fn;

%just need the first two fields to get the size of the interpolating grids
I = inputStruct;
switch lower(p.Results.type)
    case 'clean'
        [V1,V2,~,~] = ndgrid(I.elevation,I.mu0,I.SSA,I.muS);
    case 'dirty'
        [V1,V2,~,~] = ndgrid(I.elevation,I.mu0,I.SSA,I.muS,I.dust,I.soot);
end
iX1 = I.elevation;
iX2 = I.mu0;

%N-D array for output values
nBands = length(S.bands);
C = cell(nBands,1);
for k=1:nBands
    C{k} = zeros(size(V1));
end

%files with values from prepBandAlbedoLookup
switch lower(p.Results.type)
    case 'clean'
        fc = dir(fullfile(folder,'checkClean*.mat'));
        vname = 'reflArray';
    case 'dirty'
        fc = dir(fullfile(folder,'checkDirt*.mat'));
        vname = 'lapArray';
end

%look through the checkpoint files and fill in the master output file
for k=1:length(fc)
    Data = load(fullfile(fc(k).folder,fc(k).name));
    ka = Data.ka;
    kz = Data.kz;
    t = V1==iX1(ka) & V2==iX2(kz);
    assert(nnz(t)==size(Data.(vname),2),...
        'oops nnz(t)=%d, size(dataArray,2)=%d',nnz(t),size(Data.(vname),2))
    assert(nBands==size(Data.(vname),1),...
        'oops nBands=%d, size(dataArray,1)=%d',nBands,size(Data.(vname),1))
    for b=1:nBands
        C{b}(t) = Data.(vname)(b,:);
    end
end
% scale if raw values
if strcmpi(p.Results.type,'clean')
    for b=1:nBands
        % the 1.6 and 2 are hardwired, could make as arguments
        C{b} = scaleRefl(C{b},1.6,2);
    end
end

fprintf('all %d checkpoint files loaded into master output\ngenerating interpolant ...\n',...
    length(fc))

F = cell(length(C),1);
X = cell(nF,1);
for k=1:length(X)
    X{k} = I.(fn{k});
end
for b=1:length(C)
    %transform the input
    scaleX = cell(size(X));
    for k=1:length(X)
        scaleX{k} = bandScaling(S.sensor,b,fn{k},X{k});
    end
    F{b} = griddedInterpolant(scaleX,C{b},p.Results.method,p.Results.extrapolation);
end

% name the cell vectors and set the filename depending on whether clean or dirty analysis
switch type
    case 'clean'
        filename = [S.sensor '_Lookup.mat'];
        FClean = F;
        SClean = S;
        if exist(filename,'file')==2
            save(filename,'SClean','FClean','-append');
        else
            save(filename,'SClean','FClean','-v7.3');
        end
    case 'dirty'
        filename = [S.sensor '_Lookup.mat'];
        FDirty = F;
        SDirty = S;
        if exist(filename,'file')==2
            save(filename,'SDirty','FDirty','-append');
        else
            save(filename,'SDirty','FDirty','-v7.3');
        end
    otherwise
        warning('type ''%s'' not recognized so no output written',type)
end
end

function sRefl=scaleRefl(reflArray,thresholdToFix,maxValue)
% some reflArray values will be out of range when mu0 is small and muS is
% large, so fix

t = reflArray>thresholdToFix;
sRefl = reflArray;
if any(t)
    maxToFix = max(reflArray(:));
    a = thresholdToFix*(maxToFix-thresholdToFix);
    b = (maxValue-thresholdToFix)/(maxToFix-thresholdToFix);
    sRefl(t) = a+b*(sRefl(t));
end
end