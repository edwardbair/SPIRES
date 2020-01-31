function consolidate_Topography(demfile, resolution)
% consolidate_Topography(demfile, resolution)
%   Use the results from postprocess_Horizons to arrange elevations, slopes
%   and azimuths, horizons and view factors into spatial images in an HDF 5
%   file
%
% Input
%   demfile - full path to the original DEM file
%   resolution - to add to the file name, e.g. '1arcsec', or '500m'
%
% Output
%   written to AzureFileShare, or to folder argument, in HDF 5 format
%   Elevation, slope, aspect, view factor, and horizons are all in one file

if nargin<2
    if isdeployed
        disp(['usage: ' mfilename ' demfile resolution'])
    else
        disp(['usage: ' mfilename '(demfile, resolution)'])
    end
    return
end

[folder,demname,ext] = fileparts(demfile);
if isempty(folder)
    folder = '.';
end

% file names of all the filtered blocks
% input horizon folder and output diary file
f = dir(fullfile(folder,[demname '*Hfilter*.mat']));
for k=1:length(f)
    filename{k,1} = f(k).name; %#ok<AGROW>
end
disp('filtered horizon block files')
disp(filename)

% DEM
ZS = load(fullfile(folder,[demname ext]));

% output
assert(ischar(resolution),'resolution must be a character string')
outputfile = [demname '_' resolution '_Topography' '.h5'];
topofile = fullfile(folder,outputfile);

% input files
gridsize = size(ZS.Z);
RefMatrix = ZS.RefMatrix;
S = SlopeAzimuth(ZS);
% divisors for conversion to integers
if isfield(ZS,'scale')
    maxZ = double(max(abs(ZS.Z(:))))*ZS.scale;
else
    maxZ = double(max(abs(ZS.Z(:))));
end
elevdivisor = double(intmax('int16'))/(ceil(maxZ/1000)*1000);
slopedivisor = double(intmax('int16'))/90;
aspectdivisor = double(intmax('int16'))/180;
viewdivisor = double(intmax('int16'));
horzdivisor = double(intmax('int16'));

% first file to set the variables
H = load(fullfile(folder,filename{1})); % variables sinH & V are loaded
sH = H.sinH;
nHorz = size(sH,1);
viewfactor = H.V;
for k=2:length(filename)
    H = load(fullfile(folder,filename{k}));
    sH = cat(2,sH,H.sinH);
    viewfactor = cat(2,viewfactor,H.V);
end
% check that sizes match
assert(numel(viewfactor)==numel(ZS.Z),...
    'viewfactor vector and elevation grid sizes do not match');
assert(size(sH,2)==numel(ZS.Z),...
    'number of horizon sets and elevation grid size do not match');

% reshape
sH = reshape(sH',gridsize(1),gridsize(2),nHorz);
viewfactor = reshape(viewfactor',gridsize(1),gridsize(2));
V = float2integer(viewfactor,viewdivisor,0,'uint16');
horzazm = linspace(-180,180,nHorz)';

% slopes and aspects
slope = float2integer(S.Slope,slopedivisor,0,'int16');
aspect = float2integer(S.Aspect,aspectdivisor,0,'int16');
if isfield(ZS,'scale')
    Z = float2integer(double(ZS.Z)*ZS.scale,elevdivisor,0,'int16');
else
    Z = float2integer(double(ZS.Z),elevdivisor,0,'int16');
end

% write everything to an HDF 5 file
dsName = {'elevation', 'slope', 'aspect', 'viewfactor'}';
divisor = [elevdivisor slopedivisor aspectdivisor viewdivisor];
dsType = {class(Z),class(slope),class(aspect),class(V)};
deflateLevel = 9;

% remove the output file if it already exists
if exist(topofile,'file')
    delete(topofile);
end
% topography first - elevation, slope, aspect, view factor
for k=1:length(dsName)
    FillValue = intmin(char(dsType{k}));
    h5create(topofile,['/Grid/' char(dsName{k})],size(Z),...
        'Deflate',deflateLevel,...
        'ChunkSize',size(Z),...
        'DataType',char(dsType{k}),...
        'FillValue',FillValue);
    h5writeatt(topofile,['/Grid/' char(dsName{k})],'divisor',divisor(k));
end

% horizons
datasetname = '/Grid/horizons';
h5create(topofile,datasetname,size(sH),...
    'Deflate',deflateLevel,...
    'ChunkSize',[size(Z,1) size(Z,2) 1],...
    'DataType',class(sH));
h5writeatt(topofile,'/Grid/horizons','nHorizons',nHorz);
h5writeatt(topofile,'/Grid/horizons','AzmLimits',[-180 180]);
h5writeatt(topofile,datasetname,'azimuths',horzazm);
h5writeatt(topofile,'/Grid/horizons','Units','sine horizon');
h5writeatt(topofile,'/Grid/horizons','divisor',horzdivisor);

% projection attributes related to elevation dataset
switch ZS.gridtype
    case 'projected'
        h5writeProjection(topofile,'/Grid',ZS.Projection);
        h5writeatt(topofile,'/Grid','ReferencingMatrix',RefMatrix);
    case 'geographic'
        h5writeatt(topofile,'/Grid','mapprojection',ZS.gridtype);
        h5writeatt(topofile,'/Grid','ReferencingMatrix',RefMatrix);
        h5writeatt(topofile,'/Grid','geoid',ZS.Geoid);
    case 'geolocated'
        h5writeatt(topofile,'/Grid','mapprojection',ZS.gridtype);
        h5writeatt(topofile,'/Grid','geoid',ZS.Geoid);
        latitude = ZS.Lat;
        longitude = ZS.Lon;
        newds = {'/Grid/latitude' '/Grid/longitude'};
        for k=1:length(newds)
            h5create(topofile,newds{k},size(Z),...
                'Deflate',deflateLevel,...
                'ChunkSize',size(Z),...
                'DataType',class(latitude));
        end
    otherwise
        error('elevation gridtype must be ''projected'', ''geographic'', or ''geolocated''')
end

% finally, write the data
h5write(topofile,['/Grid' '/elevation'],Z);
h5write(topofile,['/Grid' '/slope'],slope);
h5write(topofile,['/Grid' '/aspect'],aspect);
h5write(topofile,['/Grid' '/viewfactor'],V);
if strcmpi(ZS.gridtype,'geolocated')
    h5write(topofile,'/Grid/latitude',latitude);
    h5write(topofile,'/Grid/longitude',longitude);
end
h5write(topofile,['/Grid' '/horizons'], sH)
disp(['consolidated topographic file ' topofile])

end