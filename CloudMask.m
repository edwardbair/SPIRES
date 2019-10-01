function out=CloudMask(R,topofile,solarAzimuth,solarZenith,sensor)
%cloud masking with Texture and SAM
%modified by NB using code from Stillinger (2019)
%input: x - Rc, surf reflectance, m x n x b

x=double(R);

% name       | MODIS              |  Landsat 8 OLI
%----------------------------------------------------
% coast. aer.|  n/a               | 1 (0.433-0.453 um)
% blueband   |  3 (0.459-0.479 um)| 2 (0.440-0.515 um)
% greenband  |  4 (0.545-0.565 um)| 3 (0.525-0.600 um)
% redband    |  1 (0.620-0.670 um)| 4 (0.630-0.680 um)
% nirband    |  2 (0.841-0.876 um)| 5 (0.845-0.885 um)
% swirband   |  6 (1.628-1.652 um)| 6 (1.560-1.660 um)
% swirband2  |  7 (2.105-2.155 um)| 7 (2.100-2.300 um)
% cirrus band|  n/a               | 9 (1.360-1.390 um)


switch sensor
    case 'MODIS'
        superpixelsize = 1e5;
        m = 0.1; 
        SAMconst = 2000;
        nItr = 5;
        blocksize=[10 10];
        bands.nirband=2;
        bands.swirband=6;
        bands.swirband2=7;
        bands.redband=1;
        bands.greenband=4;
        bands.blueband=3;
    case 'OLI'
        superpixelsize = 100;
        m = 0.5; 
        SAMconst = 200;
        nItr = 5;
        blocksize=[10 10];
        bands.nirband=5;
        %removed band 1 b/c of low bias
        bands.swirband=6;
        bands.swirband2=7;
        bands.redband=4;
        bands.greenband=3;
        bands.blueband=2;
        bands.cirrusband=9;
end

%image texture
fillpx=isnan(x);
x(fillpx)=0;
o=runGabor(0:45:135,3,1,0.5,x);
o(fillpx)=0;
o2=rssq(o,3);
T=imadjust(o2,stretchlim(o2),[]);
texthres = 0.25;
tests.notcloud = T > texthres;

%super pixels
L=SAMSuperPixels(superpixelsize,x,m,SAMconst,nItr);
validPixels = ~any(fillpx,3);
%erode by 1 to avoid texture flags along border
SE=strel('square',8); %set this value corretly
validPixels = imerode(validPixels,SE);
L(~validPixels | L < 0)=0;
tests.L = L;

x(fillpx)=-9999;
out=classifySuperPixels(tests,x,sensor,topofile,bands,...
    solarAzimuth,solarZenith,blocksize);

end


