function newR = postings2cells(R)
% newR = postings2cells(R)
%changes the raster interpretation of a raster reference object (either Map
%or Geographic) from 'postings' to 'cells'
%
%Although some data that we use come in 'postings', for example digital
%elevation models, usually 'cells' is the more appropriate raster
%interpretation in our applications.

% if already cells, just return
if strcmpi(R.RasterInterpretation,'cells')
    newR = R;
    % otherwise, change
else
    W = worldFileMatrix(R);
    if contains(class(R),'Map') && contains(class(R),'Reference')
        newR = maprasterref(W,R.RasterSize,'cells');
        if any(contains(fieldnames(R),'ProjectedCRS'))
            newR.ProjectedCRS = R.ProjectedCRS;
        end
    elseif contains(class(R),'geographic','IgnoreCase',true)
        newR = georasterref(W,R.RasterSize,'cells');
        if any(contains(fieldnames(R),'GeographicCRS'))
            newR.GeographicCRS = R.GeographicCRS;
        end
    else
        error('class(R)=%s, not recognized',class(R))
    end
end
end