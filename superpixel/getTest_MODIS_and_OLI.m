function [fullMODIS, subsetMODIS, fullOLI, subsetOLI, truth] = getTest_MODIS_and_OLI(truthFilename,MODISfilename,LS8C1_folder,runMODISFlag,modelMODISFlag)
%SNOWCLOUDERRORSTATS Calculate all the stats for CH1 of my thesis. Options
%to extract various spectral libraires from each scene (truth data or error
%class spectra
%   fullXXXX - MODIS and OLI in same projection and extent of OLI Path Row
%   aquisition
%   subsetXXX - MODIS and OLI in same projection and extent of 1000x1000
%   validation subimage from SPARCS

% MAKE OPTION TO GRAB MODIS OR COARSEN LANDSAT TO MODIS RESOLUTION, OR JUST
% USE LANDSAT

%LOAD TRUTH DATA
geotiffFilename = truthFilename;
geotiffFilename(end-7:end)='data.tif';
truthSceneInfo = geotiffinfo(geotiffFilename);
data = imread(truthFilename);
truth.snow = data == 3;
truth.cloud = data == 5;
truth.neither = data == 1 | data == 2 | data == 4;
truth.sceneInfo = truthSceneInfo;
x = truthSceneInfo.CornerCoords.X;
y = truthSceneInfo.CornerCoords.Y;

% LOAD LANDSAT 8 DATA
X = GetLandsat8_Collection1( LS8C1_folder, 'allbands' );
[LS8_15m,~,panBandInfo] = GetLandsat8_Collection1( LS8C1_folder, 'band8' ); %panband
[bqa,~,fullSceneInfo] = GetLandsat8_Collection1( LS8C1_folder, 'BQA' );
LS8_30m = X.OLI_30m;
operationalCloudMask = bqa.cloud;


% Get scene data for truth mask from full landsat scene
[Rtoa,~,~] = extractSubscene(fullSceneInfo.RefMatrix,x,y,LS8_30m,'usgs');
[panBand,~] = extractSubscene(panBandInfo.RefMatrix,x,y,LS8_15m,'usgs');
[opscloudMask,~,~] = extractSubscene(fullSceneInfo.RefMatrix,x,y,operationalCloudMask,'usgs');

fullOLI.bands = LS8_30m;
fullOLI.cloudMask = operationalCloudMask;
fullOLI.panBand =LS8_15m;

subsetOLI.bands = Rtoa;
subsetOLI.cloudMask = opscloudMask;
subsetOLI.panBand = panBand;




% can model or reproject MODIS
if runMODISFlag
    %LOAD MOD09GA DATA - cloud state at 1km / vis bands at 500m
    fileTable = infoFromMODISfilename( MODISfilename );
    info = validateHDF( MODISfilename );
    modisIn = GetMOD09GA( MODISfilename, 'allBands' );
    MOD_state1km = GetMOD09GA( MODISfilename, 'state' );
    
    %cloud - 0 clear, 1 cloudy, 2 mixed, 3 not sure
    %%%%%cloud = MOD_state1km.cloud ==1 | MOD_state1km.cloud == 2;% | MOD_state1km.cloud == 3;
    cloud = MOD_state1km.cloud;
    %cirrus - 0 none, 1 small, 2 average, 3 high
    %%%%%cirrus = MOD_state1km.cirrus == 2 | MOD_state1km.cirrus == 3;
    cirrus = MOD_state1km.cirrus;
    boundaryIdx = true(size(cloud));
    
    %REPROJECT MODIS TO UTM (working in Landsat 8 projection)
    
    %MODIS Projection Info
    
    %500m bands
    tile = char(fileTable.tile);
    [ refMatrixMODIS, InProj, InRR, ~, ~ ] = sinusoidProjMODtile( tile );
    ref500m = refMatrixMODIS.RefMatrix_500m;
    ref1km = refMatrixMODIS.RefMatrix_1km;
    R500m=ref500m;
    sz500m=size(modisIn);
    RR500m=refmatToMapRasterReference(R500m,sz500m(1:2));
    RR500m.ColumnsStartFrom='north';
    
    %1km cloud mask
    cloudIn = cat(3,cloud,cirrus,boundaryIdx);
    R1km=ref1km;
    sz1km=size(cloudIn);
    RR1km=refmatToMapRasterReference(R1km,sz1km(1:2));
    RR1km.ColumnsStartFrom='north';
    
    
    %reproject and subset to size of full landsat scene, retain modis pixel size
    
    %landsat full scene projection info
    OutProj_full = geotiff2mstruct(fullSceneInfo);
    
    [modisOut,~,~] = rasterReprojection(modisIn,InRR.RasterReference_500m,InProj,OutProj_full,...
        'XLimit',fullSceneInfo.SpatialRef.XWorldLimits,...
        'YLimit',fullSceneInfo.SpatialRef.YWorldLimits,'pixelsize',InRR.RasterReference_500m.CellExtentInWorldX);
    
    
    %[cloudOut,~,~] = rasterReprojection(cloudIn,RR1km,InProj,OutProj_full,...
    [cloudOut,~,~] = rasterReprojection(cloudIn,InRR.RasterReference_1km,InProj,OutProj_full,...
        'XLimit',fullSceneInfo.SpatialRef.XWorldLimits,...
        'YLimit',fullSceneInfo.SpatialRef.YWorldLimits,...
        'pixelsize',InRR.RasterReference_500m.CellExtentInWorldX,'method','nearest');
    
    
    %reproject and subset to size of truth mask 1000x1000 landsat subscene, retain modis pixel size
    
    %landsat subscene scene projection info
    OutProj_sub = geotiff2mstruct(truthSceneInfo);
    
    [modisOut_sub,~,~] = rasterReprojection(modisIn,InRR.RasterReference_500m,InProj,OutProj_sub,...
        'XLimit',truthSceneInfo.SpatialRef.XWorldLimits,...
        'YLimit',truthSceneInfo.SpatialRef.YWorldLimits,'pixelsize',InRR.RasterReference_500m.CellExtentInWorldX);
    
    [cloudOut_sub,~,~] = rasterReprojection(cloudIn,InRR.RasterReference_1km,InProj,OutProj_sub,...
        'XLimit',truthSceneInfo.SpatialRef.XWorldLimits,...
        'YLimit',truthSceneInfo.SpatialRef.YWorldLimits,...
        'pixelsize',InRR.RasterReference_500m.CellExtentInWorldX,'method','nearest');
    
    fullMODIS.bands = modisOut;
    fullMODIS.cloudMask = cloudOut(:,:,1);
    fullMODIS.cirrusMask = cloudOut(:,:,2);
    fullMODIS.boundary = cloudOut(:,:,3);
    
    subsetMODIS.bands = modisOut_sub;
    subsetMODIS.cloudMask = cloudOut_sub(:,:,1);
    subsetMODIS.cirrusMask = cloudOut_sub(:,:,2);
    subsetMODIS.boundary = cloudOut_sub(:,:,3);
    
elseif modelMODISFlag 
    InRR = truthSceneInfo.SpatialRef;
    %MODISfilename is just the tile... HACK
    InProj = geotiff2mstruct(truthSceneInfo);
    tile = MODISfilename;
    [ refMatrixMODIS, OutProj, OutRR, ~, ~ ] = sinusoidProjMODtile( tile );   
    [MODISmodel,~,~] = rasterReprojection(Rtoa,InRR,InProj,OutProj,...
        'pixelsize',OutRR.RasterReference_500m.CellExtentInWorldX);   
    [truth.MODIScloud,~,~] = rasterReprojection(truth.cloud,InRR,InProj,OutProj,...
        'pixelsize',OutRR.RasterReference_500m.CellExtentInWorldX,'method','nearest');  
    [truth.MODISsnow,~,~] = rasterReprojection(truth.snow,InRR,InProj,OutProj,...
        'pixelsize',OutRR.RasterReference_500m.CellExtentInWorldX,'method','nearest');  
    fullMODIS=[];
    subsetMODIS.bands = MODISmodel;   
else   
    fullMODIS=[];
    subsetMODIS=[];    
end
end

