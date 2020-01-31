function [BestPixel, SensorZenith, SolarZenith, watermask, badmask, bp] = modis8dayPixelRank(filenamesMOD09GA,cloudmask)
% [BestPixel SensorZenith SolarZenith watermask] = modis8dayPixelRank(filenamesMOD09GA,cloudmask)
%Select the best pixel from 8 MOD09GA files

% INPUT
% filenamesMOD09GA - filenames for 8 MOD09GA files
% cloudmask - 3d array with cloud mask
%
% OUTPUT
% BestPixel - 2d grid with index for day
% SensorZenith - 3d array with the sensor zenith angle
% SolarZenith - 3d array with the solar zenith angle
% watermask - 2d grid where deep ocean or continental/moderate ocean
% badmask - 3d array of bad data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Karl Rittger
% Jet Propulsion Laboratory & Earth Research Institute, UCSB
% October 8, 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate, read and scale layers from MOD09GA
sizea=[2400 2400];
cnt=length(filenamesMOD09GA);
dataset1km={'state_1km_1';'SensorZenith_1';'SolarZenith_1'};
dataset500m={'num_observations_500m';'QC_500m_1'};
datagroups={'1km';'1km';'1km';'500m';'500m'};
SensorZenith=NaN([sizea cnt],'single');
SolarZenith=NaN([sizea cnt],'single');
num_obs=NaN([sizea cnt],'single');
cloud_shadow=false([sizea cnt]);
cirrus=uint8(false([sizea cnt]));
aerosol=uint8(false([sizea cnt]));
water=uint8(false([sizea cnt]));
modland=uint8(false([sizea cnt]));
atmosphere=uint8(false([sizea cnt]));
b5qc=uint8(false([sizea cnt]));
% loop through the 8 or less files
parfor n=1:cnt
   cstruct=parload('MOD09GA_cstruct','cstruct');
   % Read
   [grid1km,grid500m] = read_MOD09GA(filenamesMOD09GA{n},...
       datagroups,[dataset1km; dataset500m]); 
   
   % Scale
   SensorZenith(:,:,n) = imresize(scale_MOD09GA(...
       grid1km.(dataset1km{2}),dataset1km{2},cstruct),sizea,'nearest');
   SolarZenith(:,:,n) = imresize(scale_MOD09GA(...
       grid1km.(dataset1km{3}),dataset1km{3},cstruct),sizea,'nearest');
   state_1km = scale_MOD09GA(grid1km.(dataset1km{1}),dataset1km{1},cstruct);
   num_obs(:,:,n) = scale_MOD09GA(grid500m.(dataset500m{1}),dataset500m{1},cstruct);
   QC_500m = scale_MOD09GA(grid500m.(dataset500m{2}),dataset500m{2},cstruct);
   
   % Layers of interest
   cloud_shadow(:,:,n)=imresize(state_1km.cloud_shadow,sizea,'nearest');
   cirrus(:,:,n)=imresize(state_1km.cirrus_detected,sizea,'nearest');
   aerosol(:,:,n)=imresize(state_1km.aerosol_quantity,sizea,'nearest');
   water(:,:,n)=imresize(state_1km.land_water_flag,sizea,'nearest');
   modland(:,:,n)=QC_500m.modland;
   atmosphere(:,:,n)=QC_500m.atmosphere;
   b5qc(:,:,n)=QC_500m.band(:,:,5);
end
% Rank pixel
bp=zeros([sizea cnt],'uint8');

% High sensor viewing angle
bp(SensorZenith>45)=2;

% High solar viewing angle
bp(SolarZenith>75)=3;

% Atmospheric correction performed
bp(atmosphere==0)=4;

% High cirrus pixels
bp(cirrus==4)=5;

% Cloud shadow
bp(cloud_shadow)=6;

% Climatology aerosol (aerosol adjustment not performed)
bp(aerosol==1)=7;

% High aerosol
bp(aerosol==4)=8;

% "Good" pixels
bp(bp==0)=9;

% Set bad data to lowest rank
% Corrected product not produced due to cloud effects all bands
bp(modland==3)=0;

% Corrected product not produced due to other reasons some or all bands may
%   be fill value
bp(modland==4)=0;

% Band 5 striping
bp(b5qc==2)=0;%Nosiy detector
bp(b5qc==3)=0;%Dead detector

% Orbital gaps
bp(num_obs==0)=0;

% Optically thick cloud pixels
bp(cloudmask)=0;

% Create a mask as well because we don't want these in output
% This should be re considered - if we select a pixel and then mask it we
% could probably just be selecting a lower ranked pixel
badmask = modland==3 | modland==4 | b5qc==2 | b5qc==3 | num_obs==0 | cloudmask;

% does the logic above actually select only some pixs?

% Find the best pixel. If there are duplicate ranks select the one with the
% smallest sensor zenith
BestPixel=zeros(sizea,'uint8');
sizea1=sizea(1);
sizea2=sizea(2);
parfor n=1:sizea1
%for n=1:sizea1
    sensz=squeeze(SensorZenith(n,:,:));%A matrix
    for m=1:sizea2
        v=squeeze(bp(n,m,:));
        bestp=max(v)==v;
        if sum(bestp)==1;
            BestPixel(n,m)=find(bestp);
        else
            sz=sensz(m,:)';
            BestPixel(n,m)=find(sz==min(sz(bestp)),1,'first');
        end
    end
end

% Land/Water flag: (could update to use GLOBCOVER to have consitent mask)
wmask=false([sizea cnt]);
wmask(water==8)=1;%deep ocean
wmask(water==7)=1;%continental/moderate ocean
wmask(num_obs==0)=0;
wmask=sum(wmask,3);
watermask=false(sizea);
watermask(wmask>=1)=1;

% Fix SensZ and SolZ
SensorZenith(num_obs==0)=NaN;
SolarZenith(num_obs==0)=NaN;
