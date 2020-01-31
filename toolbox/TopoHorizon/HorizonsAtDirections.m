function HorizonsAtDirections(demfile, varargin)
% HorizonsAtDirections(demfile [, nHorz, firstHorz, perNode] )
%   reads elevation file from file system and calculates horizons
%
%Input
%   demfile - .mat full-path filename where DEM is stored, perhaps created
%       with GridInfo, but just the fields in the structure, not the structure itself
%       (so S=load(demfile))
%       Output files go to the same folder, so you must have write
%       permission there.
%Optional input
%   nHorz - number of horizon directions, default 64, but adjusted
%       to multiples of 4 otherwise, must be specified if following arguments
%       are used, mainly for processing in a multi-node cluster)
%   firstHorz - first horizon to calculate, and its reverse, default is 1
%   perNode - number of horizons to calculate per node
%       (set this number as the increment in the parametric sweep GUI,
%       default is 1 if firstHorz is specified, all otherwise)
%
% Output
%       set of files with names of horizon files, written to same folder
%       where demfile lives (units are sines of horizon angles)
% NOTE
% The way azimuths are handled assumes that the origin of the grid is the
% NW corner. Corresponding with a right-hand coordinate system, 0 is south,
% + to east, - to west. First horizon is direction -180.

% defaults
defaultH = 64;
planet = 'earth';
minargs = 1;
maxargs = 5;
narginchk(minargs,maxargs)

if nargin<minargs
    if isdeployed
        disp(['usage: ' mfilename ' demfile [nHorz firstHorz perNode]'])
    else
        disp(['usage: ' mfilename '(demfile [, nHorz, firstHorz, perNode])'])
    end
    return
end

p = inputParser;
addRequired(p,'demfile',@ischar);
% if deployed, all args are characters
if isdeployed
    addOptional(p,'nHorz',num2str(defaultH),@ischar)
    addOptional(p,'firstHorz','',@ischar)
    addOptional(p,'perNode','',@ischar)
else
    addOptional(p,'nHorz',defaultH,@isnumeric)
    addOptional(p,'firstHorz',[],@isnumeric)
    addOptional(p,'perNode',[],@isnumeric)
end
parse(p,demfile,varargin{:});

[folder,demname,~] = fileparts(p.Results.demfile);
if isempty(folder)
    folder = '.';
end
ZS = load(p.Results.demfile);

% number of horizons
if isnumeric(p.Results.nHorz)
    N = p.Results.nHorz;
else
    N = str2double(p.Results.nHorz);
end
if isempty(N) || N==0
    N = defaultH;
elseif N>0 && N<4
    N = 4;
elseif mod(N,4)~=0
    N = 4*ceil(N/4);
end
% number of horizon pairs per node
if isempty(p.Results.perNode)
    if isempty(p.Results.firstHorz)
        iH = N/2;
    else
        iH = 1;
    end
elseif isnumeric(p.Results.perNode)
    iH = p.Results.perNode;
else
    iH = str2double(p.Results.perNode);
end
% which azimuth to start with (and go for perNodes azimuths)
if isempty(p.Results.firstHorz)
    wH = 1;
elseif isnumeric(p.Results.firstHorz)
    wH = p.Results.firstHorz;
else
    wH = str2double(p.Results.firstHorz);
end

% output diary file
if isdeployed
    diaryFolder = fullfile(folder,'output');
    if exist(diaryFolder,'file')~=7
        mkdir(diaryFolder);
    end
    StartAzureDiary(mfilename,diaryFolder,demname,...
        ['firstHorz' num2str(wH,'%04d')]);
end
disp(['file ' p.Results.demfile ' loaded'])
disp(ZS)

% angles - need just southern half of directions because horz2d computes
% forward & backward
ang = linspace(-90,90,1+N/2);
ang(1) = []; % 90 and -90 gotten on same pass

% beginning and end horizon azimuth to calculate
if wH<1
    wH = 1;
elseif wH>N/2
    wH = N/2;
end

if iH<1
    iH = 1;
elseif iH+wH-1>N/2
    iH = N/2+1-wH;
end

% other information from the input elevation structure
if isfield(ZS,'Geoid')
    Horz.Geoid = ZS.Geoid;
end

% calculate lat & lon
Z = double(ZS.Z)*ZS.scale;
switch ZS.gridtype
    case 'geographic'
        Horz.RefMatrix = ZS.RefMatrix;
        [lon, lat] = pixcenters(ZS.RefMatrix,size(Z),'makegrid');
        if isfield(ZS,'rasterReference')
            Horz.RasterReference = ZS.rasterReference;
        elseif isfield(ZS,'RasterReference')
            Horz.RasterReference = ZS.RasterReference;
        else
            Horz.RasterReference = refmatToGeoRasterReference(ZS.RefMatrix,size(ZS.Z));
        end
    case 'projected'
        Horz.RefMatrix = ZS.RefMatrix;
        Horz.Projection = ZS.Projection;
        [x,y] = pixcenters(ZS.RefMatrix,size(Z),'makegrid');
        try % some projects not supported by minvtran
            [lat,lon] = minvtran(ZS.Projection,x,y);
        catch
            [lat,lon] = projinv(ZS.Projection,x,y);
        end
        if isfield(ZS,'rasterReference')
            Horz.RasterReference = ZS.rasterReference;
        elseif isfield(ZS,'RasterReference')
            Horz.RasterReference = ZS.RasterReference;
        else
            Horz.RasterReference = refmatToMapRasterReference(ZS.RefMatrix,size(ZS.Z));
        end
    case 'geolocated'
        lat = ZS.Lat;
        lon = ZS.Lon;
    otherwise
        error('unrecognized gridtype in elevation structure')
end

% calculate horizons in the forward and backward directions for each
% azimuth
Horz.gridtype = ZS.gridtype;

for h=wH:wH+iH-1
    [sinF,sinB,azmF,azmB] =...
        Horizons2Directions(planet,ang(h),lat,lon,Z);
    sinF = medfilt2(sinF,[5 5]);
    sinB = medfilt2(sinB,[5 5]);
    % horizon files
    hF = [demname '_' 'H' num2str(h+N/2,'%04d') '.mat'];
    hB = [demname '_' 'H' num2str(h,'%04d') '.mat'];
    SaveHorz(folder,hF,hB,Horz,sinF,sinB,azmF,azmB);
    disp(['saved horizons at rotation ' num2str(ang(h)) ': files ' hF ', ' hB])
end

end

function SaveHorz(ofold,hF,hB,Horz,sinF,sinB,azmF,azmB )
% SaveHorz(ofold,hF,hB,Horz,sinF,sinB,azmF,azmB )
% writes temporary horizon information
%
% ofold - folder to write to
% hF - filename for forward horizons
% hB - filename for backward horizons
% Data - each is size of elevation grid
%   sinF - sine(forward horizons)
%   sinB - sine(backward horizons)
%   azmF - azimuths for forward horizons
%   azmB - azimuths of backward horizons

Horz.sinH = single(sinF);
Horz.azm = single(azmF);

B = whos('Horz');
if B.bytes>2^31
    save(fullfile(ofold,hF),'-struct','Horz','-v7.3');
else
    save(fullfile(ofold,hF),'-struct','Horz');
end
Horz.sinH = single(sinB);
Horz.azm = single(azmB);
B = whos('Horz');
if B.bytes>2^31
    save(fullfile(ofold,hB),'-struct','Horz','-v7.3');
else
    save(fullfile(ofold,hB),'-struct','Horz');
end

end