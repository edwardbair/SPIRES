function B = interpolateRaster(A,RR,Xq,Yq,method)
%interpolate raster and keep NaNs from propagating
%Input
%   A - input values
%   RR - input raster reference (will be either 'map' or 'geographic')
%   Xq,Yq - coordinates of interpolated values (if RR is 'geographic', then
%       Xq is longitude and Yq is latitude)
%   method - interpolation method
%Output
%   B - interpolated raster

% mapinterp or geointerp
if contains(class(RR),'Map')
    useMapInterp = true;
elseif contains(class(RR),'Geographic')
    useMapInterp = false;
else
    error('input raster reference must be ''Map'' or ''Geographic'', it''s %s',...
        clsss(RR))
end

%memory check
    RasterReprojectionMemoryCheck((numel(Xq)+numel(Yq))*8*size(A,3));

%should have enough memory, go ...
if ismatrix(A)
    if useMapInterp
        B = mapinterp(A,RR,Xq,Yq,method);
    else
        B = geointerp(A,RR,Yq,Xq,method);
    end
elseif ndims(A)==3
    for k=1:size(A,3)
        V = A(:,:,k);
        if useMapInterp
            B1 = mapinterp(V,RR,Xq,Yq,method);
        else
            B1 = geointerp(V,RR,Yq,Xq,method);
        end
        % allocate on first pass
        if k==1
            %memory check
            if ispc
                RasterReprojectionMemoryCheck(8*size(A,3)*(numel(B1)));
            end
            %allocate output space
            B = zeros(size(B1,1),size(B1,2),size(A,3));
        end
        B(:,:,k) = B1;
    end
else
    error('arrays of more than 3 dimensions not supported')
end
if all(isnan(B(:)))
    error('all interpolated values are NaN, check input and output coordinates')
end

end