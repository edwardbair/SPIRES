function [similarRGB, newcmap, newLabels] = similarRGBs(maskRGB,similarSnow,similarCloud, oldLabels, oldcmap)
%SIMILARRGBS add colors for the snow that is similar to clouds and clouds
%that are similar to snow.
%   Detailed explanation goes here

similarRGB = maskRGB;
newLabels = oldLabels;
newcmap = oldcmap;

%colors for new catagories
simSnow =  [255/255,0/255,255/255];    %3 - purple
simCloud = [0/255,200/255,0/255];   %4 - green   

%add layers to image
similarRGB(similarSnow) = 3;
similarRGB(similarCloud) = 4;

%update image info
newcmap = [newcmap;simSnow;simCloud];
newLabels(end+1)={'similarSnow'};
newLabels(end+1)={'similarCloud'};
end

