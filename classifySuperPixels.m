function out=classifySuperPixels(tests,x,sensor,topofile,bands,...
    solarAzimuth,solarZenith,blocksize)

sz=size(x);

dem=GetTopography(topofile,'elevation');
slope=GetTopography(topofile,'slope');
aspect=GetTopography(topofile,'aspect');

%put bands into variables
red=x(:,:,bands.redband);
green=x(:,:,bands.greenband);
blue=x(:,:,bands.blueband);
nir=inpaint_nans(x(:,:,bands.nirband));
swir=inpaint_nans(x(:,:,bands.swirband));
swir2=x(:,:,bands.swirband2);

sA=double(180-solarAzimuth);
sZ=double(cosd(solarZenith));

docirrus=false; %don't use cirrus band

switch sensor
    case 'OLI'
    docirrus=true;
    cirrus=x(:,:,bands.cirrusband);
    tests.cirrus= dem < 3500 & cirrus > 0.05;
    %becuase of the similar resolution of the DEM and the OLI data, I need 
    %to
    %grow the shadow mask by making the solar zeith be a little longer 
    %and the
    %solar azimuth to be plus or minus a bit.
    az_spread=5;
    z_spread=15;
    
    for i=1:2
        if i==1
            sA=double(180-solarAzimuth)-az_spread;
        elseif i==2
            sA=double(180-solarAzimuth)+az_spread;
            sZ=solarZenith+z_spread;
            sZ(sZ<0)=0;
        end
        thisshade=~GetHorizon(topofile,sA,sZ);
        terrainshade=terrainshade | thisshade;
    end
    case 'MODIS'
      terrainshade=~GetHorizon(topofile,sA,sZ);
end

NDVI=(nir-red)./(nir+red);
NDSI=(green-swir)./(green+swir);

tests.water=(NDVI < 0.01 & nir < 0.11) | ...
    (NDVI < 0.1 & NDVI > 0 & nir < 0.05);

fillednir=imfill(nir);
filledswir=imfill(swir);

deltanir=fillednir-nir;
deltaswir=filledswir-swir;

shadow_prob=min(deltanir,deltaswir);
shadows = shadow_prob > 0.02 & ~tests.water;

tests.maybecloud = terrainshade & ~shadows;
tests.maybecloudshadow = ~terrainshade & shadows; %includes shaded clouds

%haze optimized transformation - Zhang et al 2002 / Zhu and Woodcock 2012
tests.HOT=blue-0.5.*red-0.08;

%whiteness test - Gomez-Chova et al 2007 / Zhu and Woodcock 2012
vismean=(blue+green+red)./3;
tests.whiteness=(abs(blue-vismean)+abs(green-vismean)+...
    abs(red-vismean))./vismean;

%RNIRSIWR whiteness - new idea, similar to Gomez-Chova et al 2007
meansng=(swir+nir+green)./3;
whitenesssng=(abs(swir-meansng)+abs(nir-meansng)+...
    abs(green-meansng))./meansng;
tests.wsng = whitenesssng < 0.7 & blue > 0.4;

%conservative snow test - Dark soils also pass this test (thats okay) -
%expanded from Stillinger et al 2019
tests.snow = (blue > 0.6 & swir < 0.2 & swir2 < 0.2) | NDSI > 0.7;

%soil test
soilid = (swir-blue)./ (blue+swir);
soil = soilid > 0;
tests.soil = soil & blue < 0.3;  %0.25; %Validate cutoff with ecostress 
%spectral library (after PhD)

%conservative haze test - more conservative than CFMask 0.06 maybe 
%better 8/22
tests.haze = tests.HOT > 0.05;

%veg test
tests.veg = NDVI > 0.3;
tests.darkveg = NDVI > 0 & red < 0.1 ;

tests.nir=nir;

centers = regionprops('table',tests.L,'centroid');
centers = table2array(centers);

mu=sunslope(cosd(solarZenith),180-solarAzimuth,slope,aspect);

[LZAfeatures,~,~] = extractHOGFeatures(mu,...
    centers,'BlockSize',blocksize,'UseSignedOrientation',1);
frgb = cat(3,swir,nir,green);
[features,vPoints,~] = extractHOGFeatures(frgb,centers,...
    'BlockSize',blocksize,'UseSignedOrientation',1);

D = pdist2(LZAfeatures,features,'euclidean');
d = diag(D);
idx = d < 0.9;

mapPoints = round(vPoints(idx,:));

colSub=mapPoints(:,1);
rowSub=mapPoints(:,2);
ind = sub2ind([sz(1) sz(2)], rowSub, colSub);

tests.hognotcloud=false([sz(1) sz(2)]);
tests.hognotcloud(ind)=true;

%class superpixels
checksuperpixels=unique(tests.L(~tests.L==0));
boxes=regionprops(tests.L,'BoundingBox');
tests.cloudshadowL = zeros(size(tests.L));

for sp = 1:length(checksuperpixels)
    
    thisSP = checksuperpixels(sp);
    
    % Get subimage around cluster [left, top, width, height]
    BB = boxes(thisSP);
    rmin = ceil(BB.BoundingBox(2));
    rmax = rmin+BB.BoundingBox(4)-1;
    cmin = ceil(BB.BoundingBox(1));
    cmax = cmin + BB.BoundingBox(3)-1;
    
    sub.L=tests.L(rmin:rmax, cmin:cmax);
    
    if numel(sub.L) ~= 0
        frac=struct;
        idx = sub.L == thisSP;
        numPxls = nnz(idx);
        
        fn=fieldnames(tests);
        
        for ii=1:length(fn)
            sub.(fn{ii})=tests.(fn{ii})(rmin:rmax,cmin:cmax);
            frac.(fn{ii})=nnz(idx & sub.(fn{ii}))/numPxls;
        end
   
    flag.cirrus=false;
    flag.cloudfree=false;
    if docirrus
        flag.cirrus=any(idx(:) & sub.cirrus(:));
    end
        flag.hognotcloud=any(idx(:) & sub.hognotcloud(:)); 
        flag.water=median(sub.nir(:),'omitnan') < 0.1;
        flag.cloudfree=true;
        t1=frac.wsng > 0.1 && frac.notcloud < 0.1 && frac.snow < 0.1 && ...
                ~flag.water;
        t2=flag.cirrus && frac.notcloud < 0.1;
        t3=frac.haze > 0.1 && frac.notcloud < 0.1 && frac.snow < 0.1;
        if t1 || t2 || t3
            flag.cloudfree=false;
        end
    
    %remove non-cloud superpixels from the cloud mask
    if frac.water > 0.05 || frac.notcloud > 0.1 || frac.soil > 0.1 ||...
            frac.snow > 0.1 || frac.veg > 0.1 || flag.hognotcloud || ...
            frac.darkveg > 0.1 && flag.cloudfree
                sub.L(idx)=0;
        %remove possible cloud shadows
    elseif frac.maybecloudshadow >= 0.5 % 50% of super pixel is shaded 
        %where there is no terrain shade in DEM illumination
        %make sure not many bright & white pixels - could be part of a cloud
        sub.whiteness(~idx)=NaN;
        sub.red(~idx)=NaN;
        pxlwhiteness = round(sub.whiteness(idx),1);
        pixelbrightness =sub.red(idx)';
        %want lowest value whiteness and coldest temp
        pxlcharateristics = [pxlwhiteness.*-1 pixelbrightness];%THIS 
        %NEEDS TO BE IMPROVED - "Flattest reflectance pixels might be 
        %very dark- need a way tobalance illumination and brightness."
        [rows,~] = sortrows(pxlcharateristics,[1 2],'descend');
        %check 5 brightest - some rare superpixels may be very small - IF
        %there are 5 brigtht, white pixels, its not a cloud shadow, might
        %be a cloud
        [~,nP]=size(rows);
        if nP < 5
            flag.cloudshadow =  ~(all(rows(:,2) > 0.2) ...
                && all(rows(:,1) > -0.7));
        else
            flag.cloudshadow = ~(all(rows(1:5,2) > 0.2) ...
                && all(rows(1:5,1) > -0.7));
        end
        
        if flag.cloudshadow && frac.maybecloud < 0.05
            %cloud shadow
            sub.cloudshadowL(idx)=thisSP;
            sub.L(idx)=0;            
        else
            if flag.cloudfree
                sub.L(idx)=0;
            end
        end  
    end
    
    %fill full mask with superpixel classification
    tests.L(rmin:rmax, cmin:cmax) = sub.L;
    tests.cloudshadowL(rmin:rmax, cmin:cmax) = sub.cloudshadowL;
    else
            warning('zero size subimage')
    end 
end

%%%FIND CLOUD EDGES - last major issue wiht program is missing cloud edges
% find all regions around cloud boundaries an compare
%%%to cloud and other parts it touches. if more similar to cloud, add to
%%%the cloud. - IDEA FOR AFTER PHD
%simple method - grow all clouds by 1 super pixl and then shrink? - 
%IDEA FOR AFTER PHD

out.cloudmask = tests.L~=0;
out.cloudshadow = tests.cloudshadowL~=0;