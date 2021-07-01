function [ S ] = unpackLandsat8BQA( LS8_BQA, datatype )

%unpack the bit-packed Quality Assurance (QA) layer in Landsat 8
% Operatioal Land Imager (OLI) files. 
%
% replaced precollection w/ Collection 2, NB 6/8/21
%
%   from USGS Landsat 8 webstie
%   http://landsat.usgs.gov/L8QualityAssessmentBand.php
%
% Input
%    LS8_BQA - BQA variable read from GetLandsat8
%               16 bit Qa file (.tif)
%   datatype - charater sting to specify which format the landsat data is
%   in 'precollection' or 'collection1'
%
%
% Output
%   S - structure with the following fields  ( bit # for reference )
%       fill - Designated Fill QA - Bit 0
%       frame - Dropped Frame QA - Bit 1
%       terrainOcclusions - Terrain Occlusion QA - Bit 2
%       Bit 3 not populated
%       water - Water Confidence QA - Bit 4 & 5
%       cloudShadow - Reserved for Cloud Shadow QA - Bit 6 & 7
%       veg - Vegetation Confidence QA - Bit 8 & 9
%       snow_ice - Snow/Ice Confidence QA - Bit 10 & 11
%       cirrus - Cirrus Confidence QA - Bit 12 & 13
%       cloud - Cloud Confidence QA - Bit 14 & 15
%
% For the single bits (0, 1, 2, and 3):
%
%     0 = No, this condition does not exist
%     1 = Yes, this condition exists.
%
% The double bits (4-5, 6-7, 8-9, 10-11, 12-13, and 14-15) represent levels
% of confidence that a condition exists:
%
%     00 = 0 = Algorithm did not determine the status of this condition
%     01 = 1 = Algorithm has low confidence that this condition exists
%           (0-33 percent confidence)
%     10 = 2 = Algorithm has medium confidence that this condition exists
%           (34-66 percent confidence)
%     11 = 3 = Algorithm has high confidence that this condition exists
%           (67-100 percent confidence).

switch lower(datatype)
    
    case 'collection2'
        
flagname = {'fill','dilatedCloud','cirrus','cloud','cloudShadow',...
    'snow','clear','water','cloudConfidence','cloudShadowConfidence','snowIceConfidence','cirrusConfidence'};
datatype = {'logical','logical','logical','logical','logical','logical',...
    'logical','logical','uint8','uint8','uint8','uint8'};
nbits = [1 1 1 1 1 1 1 1 2 2 2 2];
bitPosition = [0 1 2 3 4 5 6 7 8 10 12 14];

    case 'collection1'
        
        flagname = {'fill','terrainOcclusion','saturation','cloud',...
    'cloudConfidence','cloudShadowConfidence','snow_iceConfidence','cirrusConfidence'};
datatype = {'logical','logical','uint8','logical','uint8','uint8',...
    'uint8','uint8'};
nbits = [1 1 2 1 2 2 2 2];
bitPosition = [0 1 2 4 5 7 9 11];

    otherwise
        error('cannot unpack LS8 data')
        
end

N = LS8_BQA; % original 16-bit integers

I = 2.^nbits-1;
for k=1:length(flagname)
        bitsToConsider = bitshift(N,-bitPosition(k));
        if strcmp(datatype{k},'logical')
            S.(flagname{k}) = bitand(bitsToConsider,I(k))>0;
        else
            S.(flagname{k}) = cast(bitand(bitsToConsider,I(k)),datatype{k});
        end
end

end