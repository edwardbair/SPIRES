function reprojecth5(h5file,target_projection,target_refmat,...
    target_rasterref,location)
% read an h5 file, reproject the data, then save it with same name appended
% with new proj name

% in: h5 file, e.g. mydata.h5
% target_projection, matlab mstruct of new projection
% target_refmat, target referencing matrix
% target_rasterref, target rastermapreference/rastergeoreference
% location, base location for datasets that you want to reproject
% in h5 file, e.g. '/Grid/MODIS_GRID_500m';

% out: h5 file with appended name
% projection info is written to the previous lower level
% ie. '/Grid' in the location above
% anything in '/' is also copied
%% 
%strip off last '/' if its there
if location(end)=='/' && length(location) > 1
    location=location(1:end-1);
end

%find location of last separator
[~,e]=regexp(location,'/');
% write projection to this location
if isempty(e)
     error('no forward slashes in your location');
elseif e(end) > 1
    base_name=location(1:e(end)-1);
elseif e(end) == 1
    base_name='/';
end
group_struct=h5info(h5file,location);
num_datasets=length(group_struct.Datasets);
%assume all datasets are the same # of row/col
sz=group_struct.Datasets(1).Dataspace.Size(1:2);

newh5file=sprintf('%s_%s.h5',h5file(1:end-3),...
    target_projection.mapprojection);

hdr = GetCoordinateInfo(h5file,group_struct.Name,sz);

%reproject all datasets in location
for i=1:num_datasets;
    datasetname=[group_struct.Name,'/',group_struct.Datasets(i).Name];
    data = h5read(h5file,datasetname);
    divisor = h5readatt(h5file,datasetname,'divisor');
    X = single(data)/divisor;
    X = reprojectRaster(X,hdr.RefMatrix,hdr.ProjectionStructure,...
        target_projection,'rasterref',target_rasterref);
    datatype=class(data);
    %eval string as function
    x=feval(datatype,X.*divisor);
    size_x=size(x);
    if length(size_x) > 2
        Chunksize=[size_x(1:2) 1];
    else
        Chunksize=size_x;
    end
    %create new dataset
    h5create(newh5file,datasetname,size(x),'Deflate',9,'ChunkSize',...
        Chunksize,'DataType',datatype)
    %write new dataset
    h5write(newh5file,datasetname,x);
    %write attributes
    write_atts(datasetname,group_struct.Datasets(i),newh5file);
end
% write projection structure to new h5 file
h5writeProjection(newh5file,base_name,target_projection);
%write Referencing Matrix
h5writeatt(newh5file,base_name,'RefMatrix',target_refmat);
%write any other attributes to '/'
base_info=h5info(h5file,'/');
if ~isempty(base_info.Attributes);
    write_atts('/',base_info,newh5file);
end
end
%%
% sub function to write attributes
function write_atts(datasetname,dataset_struct,newh5file)
num_attributes=length(dataset_struct);
for j=1:length(num_attributes);
    att_name=dataset_struct.Attributes(j).Name;
    att_value=dataset_struct.Attributes(j).Value;
    h5writeatt(newh5file,datasetname,att_name,att_value);
end
end