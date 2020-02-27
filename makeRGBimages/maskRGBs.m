function [ maskRGB, cmap, labels ] = maskRGBs( mask,varargin )
%QUICKRGBS Summary of this function goes here
% rgb matrix and cmap for a snow cloud mask.

%coloring the final mask
empty = [100/255,100/255,100/255];    %1 - grey
snow =  [102/255,194/255,165/255];    %2 - turquoiz
cloud = [255/255,217/255,47/255];   %3 - yellow   
%add no data? 
%outside = 255
 
maskRGB =uint8(zeros(size(mask.cloud)));
maskRGB(mask.neither) = 0;
maskRGB(mask.snow)=1;
maskRGB(mask.cloud)= 2;

if nargin == 2
    switch lower(varargin{1})
        case 'cirrus' 
       maskRGB(mask.cirrus) = 3;     
            
    end
    cirrus = [255/255,0/255,0/255]; %red
    cmap = [empty;snow;cloud;cirrus];
    labels = {'neither','snow','cloud','cirrus'};
else
  cmap = [empty;snow;cloud];  
  labels = {'neither','snow','cloud'};
end
end

