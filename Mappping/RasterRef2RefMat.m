function RefMat = RasterRef2RefMat(R)
%R = input raster reference object
%RefMat = output referencing matrix
[x11,y11,dx,dy,~] = cornerCoords(R);
RefMat = makerefmat(x11,y11,dx,dy);
end