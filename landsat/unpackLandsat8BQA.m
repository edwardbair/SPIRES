function [ S ] = unpackLandsat8BQA( LS8_BQA, datatype )
%[ S ] = unpackLandsat8BQA( LS8_BQA )
%unpack the bit-packed Quality Assurance (QA) layer in Landsat 8
% Operatioal Land Imager (OLI) files. 
%
% As of January 2018 this can now handle precollection or collection 1
% landsat 8 data. (https://landsat.usgs.gov/collectionqualityband)
%
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
    
    case 'precollection'
        
flagname = {'fill','frame','terrainOcclusion','reserved','water',...
    'cloudShadow','veg','snow_ice','cirrus','cloud'};
datatype = {'logical','logical','logical','logical','uint8','uint8',...
    'uint8','uint8','uint8','uint8'};
nbits = [1 1 1 1 2 2 2 2 2 2];
bitPosition = [0 1 2 3 4 6 8 10 12 14];

    case 'collection1'
        
        flagname = {'fill','terrainOcclusion','saturation','cloud',...
    'cloudConfidence','cloudShadowConfidence','snow_iceConfidence','cirrusConfidence'};
datatype = {'logical','logical','uint8','logical','uint8','uint8',...
    'uint8','uint8'};
nbits = [1 1 2 1 2 2 2 2];
bitPosition = [0 1 2 4 5 7 9 11];

    otherwise
        error('cannot unpack LS8 data - not precollection or collection 1 data')
        
end

N = LS8_BQA; % original 16-bit integers
% Timbo's code, my rewritten version is below
% for k=1:length(flagname)
%     I = 2^nbits(k)-1; % fill the lower nbits(k), leave others 0
%     if strcmp(datatype{k},'logical')
%         S.(flagname{k}) = false(size(LS8_BQA));
%         S.(flagname{k}) = bitand(N,I)>0;
%     else
%         S.(flagname{k}) = zeros(size(LS8_BQA),datatype{k});
%         S.(flagname{k}) = cast(bitand(N,I),datatype{k});
%     end
%     if k~=3
%         %move bits to read next quality value
%         N = bitshift(N,-nbits(k));
%     else
%         %skip over unassigned bit 3
%         N = bitshift(N,-2);
%     end
% end
I = 2.^nbits-1;
for k=1:length(flagname)
    %if k~=4
        bitsToConsider = bitshift(N,-bitPosition(k));
        if strcmp(datatype{k},'logical')
            %S.(flagname{k}) = bitand(N,I(k))>0;  %jeffs line - error i
            %think
            S.(flagname{k}) = bitand(bitsToConsider,I(k))>0;
        else
            S.(flagname{k}) = cast(bitand(bitsToConsider,I(k)),datatype{k});
        end
    %end
end

end