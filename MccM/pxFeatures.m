function [X,ndvi] = pxFeatures(Rs,numBands)
%PXFEATURES Feature Set for MODIS Deep Learning cloud mask
ndsi=(Rs(:,:,3)-Rs(:,:,6))./(Rs(:,:,3)+Rs(:,:,6));
ndvi=(Rs(:,:,2)-Rs(:,:,1))./(Rs(:,:,2)+Rs(:,:,1));
swirThresh=Rs(:,:,6)>0.2;
p=2;
dim=3;
pxNorm = vecnorm(Rs,p,dim)./vecnorm(ones(1,numBands)); 
%scaled 0-1 relative to max theoretical value
X=cat(3,Rs,ndsi,ndvi,swirThresh,pxNorm);
end
