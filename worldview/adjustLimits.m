function [xnew,ynew] = adjustLimits(xlimit,ylimit,pixelsize)
% pixelsize is [height width]
xnew = xlimit;
ynew = ylimit;

if mod(xlimit(1),pixelsize(2))~=0
    xnew(1) = pixelsize(2)*floor(xlimit(1)/pixelsize(2));
end
if mod(xlimit(2),pixelsize(2))~=0
    xnew(2) = pixelsize(2)*ceil(xlimit(2)/pixelsize(2));
end
if mod(ylimit(1),pixelsize(1))~=0
    ynew(1) = pixelsize(1)*floor(ylimit(1)/pixelsize(1));
end
if mod(ylimit(2),pixelsize(1))~=0
    ynew(2) = pixelsize(1)*ceil(ylimit(2)/pixelsize(1));
end
% preserve aspect ratio
% aspectRatio = (max(xlimit)-min(xlimit))/(max(ylimit)-min(ylimit));
% newAspectRatio = (max(xnew)-min(xnew))/(max(ynew)-min(ynew));
% thresh = 1.e-6;
% if abs(aspectRatio-newAspectRatio)>thresh
%     if newAspectRatio>aspectRatio
%         ynew(2) = ynew(2)+pixelsize(1);
%     else
%         xnew(2) = xnew(2)+pixelsize(2);
%     end
% end
end