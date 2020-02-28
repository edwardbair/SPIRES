function [spFeaturesTbl] = spFeatureSets(L,pixelR,bandnames,varagin)
%SPFEATURESETS Summary of this function goes here
%   Detailed explanation goes here
%inputs
%L - label matrix for which superpixel each pixel belongs to. Valid lables
%are integers grater than 0. negative values not allowed. 0 is assumed
%masked data not to be calculated.
%
%pixelR - MxNxB martix of pixel reflectance values. bands can be any number
%or combination of bands
%bandnames - cell array of the band names, will become the column headings
%in the output table
% varagin
% median values from all component pixels for the features in each input as name value pairs.
%for example ... 'canopyCover',X) X is a matrix of size MxNx1 of canopy
%cover values and the output table of median values for each superpixel
%will get the column name 'canopyCover'
%output
%spFeaturesTbl - table of output features for each superpixel
%superpixel lable from L | number of pixels in superpixel | median values
%from each band in pixelR | median values from varagin features

%find number of superpixels
checkSP = unique(L(L>0));
numSps = length(checkSP);

%0 is background(no properties calculated
invalidP=L<0;
if any(invalidP)
    warning('some invalid pixel labels... negative lables not allowed')
    L(invalidP)=0;
end


%find number of features

%inialize output
spFeatures = zeros(numSPs,numFeatures);

%update cluster centers - make its own funciton????

pixPerSP = regionprops(L,'area');

spFeatures(1,:) = 1:max(L(:));
spFeatures(2,:) =cell2mat({pixPerSP.Area});


%make work for median values later -defualt funciton is mean

%band reflectance sp means
for j = 3:bands+2
    sB = regionprops(L,X(:,:,j),'MeanIntensity');
    spFeatures(j,:) = cell2mat({sB.MeanIntensity});
end

%additonal features sp means
if numargs>3
    for j = 
        
        F = varagin();
        sB = regionprops(L,F(:,:,j),'MeanIntensity');
        spFeatures(j,:) = cell2mat({sB.MeanIntensity});
    end
end

spFeaturesTbl = array2table(spFeatures,'VariableNames',varCell);
end

