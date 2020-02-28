function [ fRGB ] = landsat8fRGB( Rs, RGB )
%LANDSAT8FRGB Summary of this function goes here
%    color matrix for plotting a browse image of a scene 
% in a figure
% rs - reflectance data for plotting
% RGB - [redBand greenBand blueBand]

%RGB = [6,5,4]; %false color
%RGB = [4,3,2];%true color
% 2,3,4 
% 4,3,2 
%flase color = 3,4,5

% Scale the surface reflectance
fRGB = NaN([size(Rs(:,:,1)) 3],'single');


for n=1:3
    tmp = single(Rs(:,:,RGB(n)));
    tmp(tmp==-9999)=NaN;
    %tmp(Landsat.cfmask==255)=0;
    fRGB(:,:,n)=tmp/10000;
end
lims=[0.01 0.99];
LOW_HIGH = stretchlim(fRGB,lims);
fRGB=imadjust(fRGB,LOW_HIGH);
fRGB=uint8(fRGB*255);
end

