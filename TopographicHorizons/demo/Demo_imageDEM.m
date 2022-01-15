function et = Demo_imageDEM(Z,R)
%images of DEM similar to Fig. 1

% intrinsic coordinates
[xI,yI] = meshgrid(1:R.RasterSize(2),1:R.RasterSize(1));
[lat,lon] = intrinsicToGeographic(R,xI,yI);

tic; % start the timer

% image of the elevations
figure('Name','Fig. 1a Elevation')
ax = setAxes(R,true); %#ok<NASGU>
geoshow(double(Z),R,'DisplayType','surface');
axis equal tight
colormap(demcmap(double(Z)))
colorbar('SouthOutside')

% shaded relief
figure('Name','Fig. 1b Shaded Relief')
ax = setAxes(R,true); %#ok<NASGU>
surflsrm(lat,lon,double(Z),[45 315])
axis equal tight

et = toc;

fprintf('This code %s reproduces the left 2 images in Fig. 1\n',...
    mfilename)

end