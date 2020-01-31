function [AVNG,varargout] = parse_AVNG_ENVI(imgFile,varargin)
% [AVNG] = parse_AVNG_ENVI(imgFile [,hdrFile])
% [AVNG,matfilename] = parse_AVNG_ENVI(imgFile [,hdrFile])
%reads AVIRIS-NG file in ENVI format, rotates so that top of image is
%toward the north, approximately, if the original rotation is >45 and <=90
%degrees
%
%Input
%   imgFile - AVIRIS-NG file in ENVI format
%Optional input
%   hdrFile - ENVI header (default is to add '.hdr' to image file name)
%
%Output
%   AVNG - structure with AVIRIS-NG image or associated variables, in BIP
%       format if radiance or reflectance image, BSQ otherwise
%Optional output
%   matfilename - writes result to file in same folder as input file(s)
%       eliminating name of structure itself, but just the variables
%

narginchk(1,2)
nargoutchk(1,2)

p = inputParser;
addRequired(p,'imgFile',@ischar)
addOptional(p,'hdrFile','',@ischar)
parse(p,imgFile,varargin{:})

% datetime (not pixel specific, see the obs or obs_ort files for that)
[~,fn,~]=fileparts(p.Results.imgFile);
allk = strfind(fn,'2'); % starts with the first 2
k = allk(1);
datepart = fn(k:end);
allk = strfind(datepart,'_'); % ends with a '_';
kend = allk(1)-1;
dt = datetime(datepart(1:kend),'InputFormat','uuuuMMdd''t''HHmmss');

% header information
if isempty(p.Results.hdrFile)
    [H,I] = read_AVNG_ENVIhdr([p.Results.imgFile '.hdr']);
else
    [H,I] = read_AVNG_ENVIhdr(p.Results.hdrFile);
end

%read input file, goes into BSQ format
X = multibandread(p.Results.imgFile,[I.lines I.samples I.bands],H.data_type,...
    I.header_offset,I.interleave,'n');
if isscalar(H.badNegative) && H.badNegative
    X(X<0) = H.ignoreValue;
elseif isfield(H,'badValue')
    X(X==H.badValue) = H.ignoreValue;
end
if ~isscalar(H.badNegative)
    for k=1:size(X,3)
        if H.badNegative(k)
            holdX = squeeze(X(:,:,k));
            holdX(holdX<0) = H.ignoreValue;
            X(:,:,k) = holdX;
        end
    end
end

%re-orient depending on rotation
if H.doMap
    if H.mapping.rotation==90
        % referencing matrix, rotated 90 deg counterclockwise, so x and y
        % reversed (x along rows, y along columns)
        origR = makerefmat(H.mapping.cornerY,H.mapping.cornerX,...
            H.mapping.pixelSize(2),H.mapping.pixelSize(1));
        % rotated coordinates
        [y,x] = pixcenters(origR,size(X,1),size(X,2),'makegrid');
        X = rot90(X);
        x = rot90(x);
        y = rot90(y);
        D = [x(2,1)-x(1,1) y(2,1)-y(1,1); x(1,2)-x(1,1) y(1,2)-y(1,1)];
        RefMatrix = makerefmat(x(1,1),y(1,1),D(:,1)',D(:,2)');
    elseif H.mapping.rotation<=45 && H.mapping.rotation>0
        RefMatrix = makerefmat(H.mapping.cornerX,H.mapping.cornerY,...
            H.mapping.pixelSize(1)*[sind(H.mapping.rotation) cosd(H.mapping.rotation)],...
            H.mapping.pixelSize(2)*[-cosd(H.mapping.rotation) sind(H.mapping.rotation)]);
    elseif H.mapping.rotation>45 && H.mapping.rotation<90
        origR = makerefmat(H.mapping.cornerX,H.mapping.cornerY,...
            H.mapping.pixelSize(1)*[sind(H.mapping.rotation) cosd(H.mapping.rotation)],...
            H.mapping.pixelSize(2)*[-cosd(H.mapping.rotation) sind(H.mapping.rotation)]);
        [x,y] = pixcenters(origR,size(X,1),size(X,2),'makegrid');
        X = rot90(X);
        x = rot90(x);
        y = rot90(y);
        D = [x(2,1)-x(1,1) y(2,1)-y(1,1); x(1,2)-x(1,1) y(1,2)-y(1,1)];
        RefMatrix = makerefmat(x(1,1),y(1,1),D(:,1)',D(:,2)');
    else % no rotation
        RefMatrix = makerefmat(H.mapping.cornerX,H.mapping.cornerY,...
            H.mapping.pixelSize(1),-H.mapping.pixelSize(2));
    end
end

% output structure
AVNG.description = H.description;
AVNG.filename = fn;
AVNG.acquisitionDateTime = dt;
if H.doMap
    AVNG.gridtype = 'projected';
    AVNG.interleave = H.outerLeave;
    AVNG.RefMatrix = RefMatrix;
    AVNG.projection = H.projection;
end
if isfield(H,'wavelength')
    AVNG.wavelength = H.wavelength;
    AVNG.wavelength_units = 'nm';
end
if isfield(H,'bandName')
    AVNG.bandName = H.bandName;
end
AVNG.data_type = H.data_type;
AVNG.ignoreValue = H.ignoreValue;
if strcmpi(AVNG.interleave,'bip')
    X = bsq2bip(X);
end

% choose name of variable based on its type
AVNG.(H.outputVariableName) = cast(X,H.data_type);

% write to file if 2nd argument
if nargout>1
    % add '.mat' to filename
    matfilename = [imgFile '.mat'];
    save(matfilename,'-struct','AVNG','-v7.3');
    varargout{1} = matfilename;
end

end