function B = interpolateRaster(geolocated,X,Y,A,Xq,Yq,method)
%interpolate raster and keep NaNs from propagating
%Input
%   geolocated - logical, true if inputs are geolocated
%   X,Y - coordinates of input values
%   A - input values
%   Xq,Yq - coordinates of interpolated values
%   method - interpolation method
%Output
%   B - interpolated raster

if ismatrix(A)
    if geolocated
        F = scatteredInterpolant(X(:),Y(:),A(:),method,'none');
        B = F(Xq,Yq);
    else
        B = interp2(X,Y,A,Xq,Yq,method,NaN);
    end
elseif ndims(A)==3
    for k=1:size(A,3)
        V = A(:,:,k);
        if geolocated
            F = scatteredInterpolant(X(:),Y(:),V(:),method,'none');
            B1 = F(Xq,Yq);
        else
            B1 = interp2(X,Y,V,Xq,Yq,method,NaN);
        end
        % allocate on first pass
        if k==1
            B = zeros(size(B1,1),size(B1,2),size(A,3));
        end
        B(:,:,k) = B1;
    end
else
    error('arrays of more than 3 dimensions not supported')
end

end

