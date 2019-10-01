function out=CloudMaskTesting(x)
%cloud masking with Texture and SAM
%Stillinger et al
x=double(x);

%set variables here
superPixelSize = 1e6;
m = 0.1; 
SAMconst = 200;
nItr = 5;

nirband=2;
swirband=6;
swirband2=7;
redband=1;
greenband=4;
blueband=3;

%image texture
fillpx=isnan(x);
x(fillpx)=0;
o=runGabor(0:45:135,3,1,0.5,x);
o(fillpx)=0;
o2=rssq(o,3);
T=imadjust(o2,stretchlim(o2),[]);
texThres = 0.25;
tests.notCloud = T>texThres;

%super pixels
sz=size(x);
L = SAMSuperPixels(superPixelSize,x,m,SAMconst,nItr);

%put bands into variables
red=Rc(:,:,redband);
green=Rc(:,:,greenband);
blue=Rc(:,:,bluenand);
nir=inpaint_nans(x(:,:,nirband));
swir=inpaint_nans(x(:,:,swirband));
swir2=x(:,:,swir3band);

NDVI=(x(:,:,nir)-x(:,:,red))./(nir+red);
NDSI=(green-swir)./(green+swir);

tests.waterTest=(NDVI < 0.01 & nir < 0.11) | ...
    (NDVI < 0.1 & NDVI > 0 & nir < 0.05);

fillednir=imfill(nir);
filledswir=imfill(swir);

deltanir=fillednir-nir;
deltaswir=filledswir-swir;

shadow_prob=min(deltanir,deltaswir);
shadows = shadow_prob > 0.02 & ~waterTest;

%haze optimized transformation - Zhang et al 2002 / Zhu and Woodcock 2012
tests.HOT=blue-0.5*red-0.08;

%whiteness test - Gomez-Chova et al 2007 / Zhu and Woodcock 2012
visMean=(blue+green+red)./3;
tests.whiteness=(abs(blue-visMean)+abs(green-visMean)+...
    abs(red-visMean))./visMean;

%RNIRSIWR whiteness - new idea, similar to Gomez-Chova et al 2007
meansng=(swir+nir+green)./3;
whitenesssng=(abs(swir-meansng)+abs(nir-meansng)+...
    abs(green-meansng))./meansng;
tests.wsng = whitenesssng<.7 & X(:,:,2)>.4;

%conservative snow test - Dark soils also pass this test (thats okay) -
%expanded from Stillinger et al 2019
tests.snow = (blue>0.6 & swir<0.2 & swir2<0.2) | NDSI>0.7;

%soil test
%soilIdxB1 = (X(:,:,swirband)-X(:,:,1))./(X(:,:,1)+X(:,:,swirband));
soilId = (swir-blue)./ (blue+swir);
soil = soilId > 0;
tests.soil = soil & blue < 0.3;  %0.25; %Validate cutoff with ecostress spectral library (after PhD)

%conservative haze test - more conservative than CFMask 0.06 maybe better 8/22
tests.hazeTest = HOT > 0.05;

%veg test
tests.veg = NDVI > 0.3;
tests.darkVeg = NDVI > 0 & red < 0.1 ;

terrainShade=GetHorizon(topofile,double(180-solarAzimuth),...
    double(solarZenith));

maybeCloud = terrainShade & ~shadows;
maybeCloudShadow = ~terrainShade&shadows; %includes shaded clouds
mu=sunslope(cosd(solarZenith),180-solarAzimuth,...
    GetTopography(topofile,'Slope'),GetTopography(topofile,'Aspect'));

centers = regionprops('table',L,'Centroid');
centers = table2array(centers);

[LZAfeatures,~,~] = extractHOGFeatures(mu,...
    centers,'UseSignedOrientation',1);
rgb = cat(3,red,green,blue);
[features,vPoints,~] = extractHOGFeatures(rgb,centers,...
    'UseSignedOrientation',1);

D = pdist2(LZAfeatures,features,'euclidean');
d = diag(D);
idx = d < 0.9;

mapPoints = round(vPoints(idx,:));

colSub=mapPoints(:,1);
rowSub=mapPoints(:,2);
ind = sub2ind([sz(1) sz(2)], rowSub, colSub);

hogNotCloud = false([sz(1) sz(2)]);
hogNotCloud(ind)=true;

%class superpixels
checkSuperPixels = unique(L(~L==0));
boxes=regionprops(L,'BoundingBox');
cloudShadowL = zeros(size(L));

for sp = 1:length(checkSuperPixels)
    
    thisSP = checkSuperPixels(sp);
    
    % Get subimage around cluster [left, top, width, height]
    BB = boxes(thisSP);
    rmin = ceil(BB.BoundingBox(2));
    rmax = rmin+BB.BoundingBox(4)-1;
    cmin = ceil(BB.BoundingBox(1));
    cmax = cmin + BB.BoundingBox(3)-1;
    
    subL = L(rmin:rmax, cmin:cmax);
    subShadowL = cloudShadowL(rmin:rmax, cmin:cmax);
    
    if numel(subL) ~= 0

        subSP = subL == thisSP;
        numPxls = nnz(subSP);
    
        sub=struct;
        fract=struct;
        fn=fieldnames(tests);
        for ii=1:length(fn)
            sub.(fn)=tests.(fn{ii})(rmin:rmax);
            fract.(fn)=nnz(subSP & sub.(fn))/numPxls;
        end
        flagHOG =  any(subSP(:) & subHOG(:))
    
    %extract superpixel spectral test results
%     subNC = notCloud(rmin:rmax, cmin:cmax);%,
%     subSnow = snow(rmin:rmax, cmin:cmax);%, 
%     subSoil = soil(rmin:rmax, cmin:cmax);
%     subWater = waterTest(rmin:rmax, cmin:cmax);
%     subCloudShadow = maybeCloudShadow(rmin:rmax, cmin:cmax);
%     subMaybeCloud = maybeCloud(rmin:rmax, cmin:cmax);
%     subwsng = wsng(rmin:rmax, cmin:cmax);
%     subHOG = hogNotCloud(rmin:rmax, cmin:cmax);
%     subVeg = veg(rmin:rmax, cmin:cmax);
%     subCirrus = cirrus(rmin:rmax, cmin:cmax);
%     subHaze = hazeTest(rmin:rmax, cmin:cmax);
%     subB4 = b4(rmin:rmax, cmin:cmax);
%     subB5 = b5(rmin:rmax, cmin:cmax);
%     subWhite = whiteness(rmin:rmax, cmin:cmax);
%     subDarkVeg = darkVeg(rmin:rmax, cmin:cmax);
    
    % fractional cover tests
%     fracCloudShadow = nnz(subSP & subCloudShadow)/numPxls;
%     fracWater = nnz(subWater&subSP)/numPxls;
%     fracMaybeCloud = nnz(subMaybeCloud&subSP)/numPxls;
%     fracNC = nnz(subSP & subNC)/numPxls; %too much texture?
%     fracSnow = nnz(subSP & subSnow)/numPxls; %def snow?
%     fracSoil = nnz(subSP & subSoil)/numPxls; %def soil?
%     fracW653 = nnz(subSP & subw653)/numPxls; %white clouds in 653 fRGB
%     fracVeg= nnz(subSP & subVeg)/numPxls; %def veg?
%     fracHaze = nnz(subSP & subHaze)/numPxls;
%     fracDarkVeg = nnz(subSP & subDarkVeg)/numPxls;
    
    %flags
%     flagCirrus=any(subSP(:) & subCirrus(:)); %def cirrus?
%     flagHOG =  any(subSP(:) & subHOG(:));%OLI texture match DEM illumination texture?
%     
    subB5(~subSP)=nan;
    waterFlag = median(subB5(:),'omitnan')<0.1;
    
    
    cloudFreeFlag = false;
    if fracW653>0.1 && fracNC<.1 && fracSnow<0.1 && ~waterFlag
        %do nothing, still cloud
    elseif flagCirrus && fracNC<.1
        %do nothing, still cloud
    elseif fracHaze>0.1 && fracNC<.1 && fracSnow<0.1
        %do nothing, still cloud
    else
        %flag as possibly not a cloud
        cloudFreeFlag = true;
    end
    
    %remove non cloud superpixels from the cloud mask
    if fracWater>0.05 || fracNC>0.1 || fracSoil>0.1 || fracSnow>0.1 || fracVeg>0.1 || flagHOG || fracDarkVeg >.1
        
        if cloudFreeFlag
            subL(subSP) =0;
        end
        
        %remove possible cloud shadows
    elseif fracCloudShadow >= 0.5 % 50% of super pixel is shaded where there is no terrain shade in DEM illumination
        
        %make sure not many bright & white pixles - could be part of a cloud
        subWhite(~subSP)=nan;
        subB4(~subSP)=nan;
        pxlWhiteness = round(subWhite(subSP),1);
        pixelBrightness =subB4(subSP);
        %want lowest value whiteness and coldest temp
        pxlCharateristics = [pxlWhiteness.*-1 pixelBrightness];%THIS NEEDS TO BE IMPROVED - "Flattest reflectance pixels might be very dark- need a way tobalance illumination and brightness."
        [R,~] = sortrows(pxlCharateristics,[1 2],'descend');
        
        %check 5 brightest - some rare superpixels may be very small - IF
        %there are 5 brigtht, white pixels, its not a cloud shadow, might
        %be a cloud
        [~,nP]=size(R);
        if nP<5
            cloudShadowFlag =  ~(all(R(:,2) > 0.2) && all(R(:,1) > -0.7));%all(R(:,1) < -0.7) && all(R(:,2) < 0.2); %
        else
            cloudShadowFlag = ~(all(R(1:5,2) > 0.2) && all(R(1:5,1) > -0.7));%all(R(1:5,1) < -0.7) && all(R(1:5,2) < 0.2);
        end
        
        if cloudShadowFlag && fracMaybeCloud<0.05
            %cloud shadow
            subShadowL(subSP) = thisSP;
            subL(subSP) =0;            
        else
            if cloudFreeFlag
                subL(subSP) =0;
            end
        end
        
    end
    
    %fill full mask with superpixel classification
    L(rmin:rmax, cmin:cmax) = subL;
    cloudShadowL(rmin:rmax, cmin:cmax) = subShadowL;
    else
            warning('0 size subimage')
    end 
end



%%%FIND CLOUD EDGES - last major issue wiht program is missing cloud edges
% find all regions around cloud boundaries an compare
%%%to cloud and other parts it touches. if more similar to cloud, add to
%%%the cloud. - IDEA FOR AFTER PHD
%simple method - grow all clouds by 1 super pixl and then shrink? - IDEA FOR AFTER PHD


cloudMask = L~=0;
cloudShadowMask = cloudShadowL~=0;
end


