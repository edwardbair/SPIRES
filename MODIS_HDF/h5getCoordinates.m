function [ hdr ] = h5getCoordinates( h5file )
% [ hdr ] = h5getCoordinates( h5file )
%

% map projection or geoid
baseloc='/Grid/';
Info = h5info(h5file,baseloc);
for k=1:length(Info.Datasets)
    Attributes = Info.Datasets(k).Attributes;
    switch Info.Datasets(k).Name
        case 'MapProjection'
            mstruct = defaultm;
            for n=1:length(Attributes)
                X = h5readatt(h5file,[baseloc 'MapProjection'],Attributes(n).Name);
                if iscolumn(X)
                    mstruct.(Attributes(n).Name) = X.';
                else
                    mstruct.(Attributes(n).Name) = X;
                end
            end
            hdr.Projection = defaultm(mstruct);
        case 'SpatialRef'
            S = struct;
            for n=1:length(Attributes)
                X = h5readatt(h5file,[baseloc 'SpatialRef'],Attributes(n).Name);
                if iscolumn(X)
                    S.(Attributes(n).Name) = X.';
                else
                    S.(Attributes(n).Name) = X;
                end
            end
            if contains(S.RasterClass,'MapCells')
                hdr.rasterref =...
                    refmatToMapRasterReference(S.RefMatrix,S.RasterSize,'cells');
            elseif contains(S.RasterClass,'MapPostings')
                hdr.rasterref =...
                    refmatToMapRasterReference(S.RefMatrix,S.RasterSize,'postings');
            elseif contains(S.RasterClass,'GeographicCells')
                hdr.rasterref = ...
                    refmatToGeoRasterReference(S.RefMatrix,S.RasterSize,'cells');
            elseif contains (S.RasterClass,'GeographicPostings')
                hdr.rasterref = ...
                    refmatToGeoRasterReference(S.RefMatrix,S.RasterSize,'postings');
            else
                error('RasterClass %s not recognized',S.RasterClass)
            end
        case 'geoid'
            for n=1:length(Attributes)
                S = struct;
                X = h5readatt(h5file,[baseloc 'SpatialRef'],Attributes(n).Name);
                if iscolumn(X)
                    S.(Attributes(n).Name) = X.';
                else
                    S.(Attributes(n).Name) = X;
                end
            end
            hdr.geoid = S.geoid;
    end
end
end