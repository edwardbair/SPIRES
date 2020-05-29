function [ hdr ] = GetCoordinateInfo( h5file,location,rastersize )
% [ hdr ] = GetCoordinateInfo( h5file,location,rastersize )
%

% 'location' is problematic, because sometimes the information needed is in
% '/Grid', but otherwise in '/Grid/500m' or '/Grid/MODIS_GRID_500m'.
% Logically, information about the projection should be in '/Grid' because
% it's independent of pixel size, but the referencing matrix has to be in
% the dataset with the right pixel size.
% So set an alternate location.
location2 = '/Grid';

rastersize = rastersize(1:2);
try
proj = h5readatt(h5file,location,'mapprojection');
projLocation = location;
catch
    proj = h5readatt(h5file,location2,'mapprojection');
    projLocation = location2;
end
if iscell(proj)
   proj=proj{:}; %problem with char stored as cell in PC
end

switch proj
    case 'geographic'
        hdr.gridtype = 'geographic';
        hdr.RefMatrix = h5readatt(h5file,location,'ReferencingMatrix');
        hdr.Geoid = h5readatt(h5file,projLocation,'geoid');
        hdr.RasterReference = refmatToGeoRasterReference(hdr.RefMatrix,rastersize);
    case 'geolocated'
        hdr.gridtype = 'geolocated';
        hdr.Geoid = h5readatt(h5file,projLocation,'geoid');
        hdr.Lat = h5read(h5file,[location '/latitude']);
        hdr.Lon = h5read(h5file,[location '/longitude']);
    otherwise
        info=h5info(h5file,projLocation);
        nFields = size(info.Attributes);
        hdr.gridtype = 'projected';
        hdr.RefMatrix = h5readatt(h5file,location,'ReferencingMatrix');
        hdr.RasterReference = refmatToMapRasterReference(hdr.RefMatrix,rastersize);
        % create projection structure
        mstruct = defaultm(proj);
        for k=1:nFields
            thisfield = info.Attributes(k).Name;
            
            if ~(strcmpi(thisfield,'mapprojection') ||...
                    strcmpi(thisfield,'referencingmatrix'))
                v=info.Attributes(k).Value;
                if iscell(v) % problem with char stored as cell in PC
                    v=v{:};
                end
                mstruct.(thisfield) = v;
            end
        end
        fn = fieldnames(mstruct);
        for k=1:length(fn)
            if iscolumn(mstruct.(fn{k}))
                mstruct.(fn{k}) = mstruct.(fn{k})';
            end
        end
        hdr.ProjectionStructure = defaultm(mstruct);
        % may need to add more here as we work with projected data
end
end