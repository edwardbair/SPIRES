function writeh5stcubes(filename,dStruct,hdr,matdate,member,Value)
% write out spacetime h5 cubes for each endmemeber,
% based on Jeff Dozier's cube2file
% %input:
%filename - filename to write out, .h5 or .mat
% dStruct - struct w/ fields:
%   divisor - to convert  using offset + value/divisor
%   dataType - numeric type, e.g. 'uint16', 'int16', 'uint8', or 'int8'
%   maxVal - maximum value
%   FillValue (optional) - null value, scalar
%   units (optional) - units for Value, string
% hdr - header w/ geographic info, see GetCoordinateInfo
% matdate - MATLAB datenums
% member - member name, string i.e. 'fsca', 'grainradius', or 'dust'
% Value - values of endmemeber

persistent already

if exist(filename,'file')==0
    X.(member)=[];
    already=[];
end

% convert to scaled integer
D = dStruct;

X.(member) = float2integer(Value,D.(member).divisor,0,...
    D.(member).dataType,0,D.(member).maxVal);
if isfield(D.(member),'FillValue')
    if D.(member).FillValue~=0
        X.(member)(isnan(Value)) = D.(member).FillValue;
    end
end

% .h5 file

group = '/Grid/MODIS_GRID_500m'; % MODSCAG data at 500 m resolution
arraySize = size(Value);
chunkSize = [arraySize(1) arraySize(2) 1];
deflateLevel = 9;
% write data to file
if isfield(D.(member),'FillValue')
    h5create(filename,[group '/' member],arraySize,...
        'Deflate',deflateLevel,...
        'ChunkSize',chunkSize,...
        'DataType',D.(member).dataType,...
        'FillValue',D.(member).FillValue)
else
    h5create(filename,[group '/' member],arraySize,...
        'Deflate',deflateLevel,...
        'ChunkSize',chunkSize,...
        'DataType',D.(member).dataType)
end
h5write(filename,[group '/' member],X.(member))
h5writeatt(filename,[group '/' member],'divisor',D.(member).divisor)
if isfield(D.(member),'units')
    h5writeatt(filename,[group '/' member],'units',D.(member).units)
end

% initial values
if isempty(already)
    already = true;
    % referencing matrix and projection information
    h5writeProjection(filename,'/Grid',hdr.ProjectionStructure)
    h5writeatt(filename,group,'ReferencingMatrix',hdr.RefMatrix)
    % write dates
    ISOdates = datenum2iso(matdate,7);
    MATLABdates = matdate;
    h5writeatt(filename,'/','MATLABdates',MATLABdates);
    h5writeatt(filename,'/','ISOdates',ISOdates);
end
end