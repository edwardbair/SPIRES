function [ maskRGB, cmap, labels ] = errorMaskRGBs( Image )
%QUICKRGBS Summary of this function goes here
% rgb matrix and cmap for error analysis results of cloud mask.


%PRINT CFMASK RESULTS AS FOLLOWING:
%FP Snow -red
%FP Other - yellow
%TP Cloud - white
%TN - black
%FN snow - cyan
%FN other - blue

%8/20/2019- changed so that TP color is not the same as initalizing zeros
%of empty matirx

%coloring the final mask with 8 bit colors
truePos = [1,1,1]; %white
falsePosSnow = [255/255,0,0]; %red
falsePosOther = [255/255,217/255,47/255]; % yellow  %[255/255,128/255,0/255]; %orange
trueNegSnow =  [102/255,194/255,165/255]; %turquoiz
trueNegOther = [0,0,0]; %black
falseNeg = [0,0,255/255]; %blue
 
%seperate snow FP and TN pixels from the rest
idxSnow = strcmp('snow',Image.cloud.key_imageTNFP);
[rS,~] =find(idxSnow);
snowNum = Image.cloud.key_imageTNFP{rS,1};

idxNeither = strcmp('neither',Image.cloud.key_imageTNFP);
[rN,~] =find(idxNeither);
otherNum = Image.cloud.key_imageTNFP{rN,1};
 

%FP just snow
snowFP = Image.cloud.imageFP == snowNum;
otherFP = Image.cloud.imageFP == otherNum;
 
%TN just snow
snowTN = Image.cloud.imageTN == snowNum;
otherTN = Image.cloud.imageTN == otherNum;

maskRGB =uint8(zeros(size(snowFP)));
maskRGB(Image.cloud.TP) = 4;
maskRGB(snowFP) =1;
maskRGB(otherFP) = 2;
maskRGB(snowTN) = 3;
maskRGB(otherTN) = 0;
maskRGB(Image.cloud.FN) = 5;




  cmap = [trueNegOther;falsePosSnow;falsePosOther;trueNegSnow;truePos;falseNeg];  
  labels = {'TN Neither','FP Snow','FP Neither','TN Snow',...
      'TP Cloud','FN Cloud'};
  
end

