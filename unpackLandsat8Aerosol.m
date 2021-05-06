function S = unpackLandsat8Aerosol( LS8_AERO )

% unpack Landsat8 aerosol
%
% Input
%    LS8_AERO - SR_QA_AEROSOL values already read from geotiff 
%
% Output
%   S - structure with the following fields  ( bit # for reference )
%       fill - Designated Fill QA - Bit 0
%       validAerosol - aerosol retrieval is valid - Bit 1
%       water - pixel is water - Bit 2
%       unused - Bit 3
%       unused - Bit 4
%       interpolatedAerosol - aerosol value is interpolated - Bit 5
%       aerosolLevel - 4 values - Bit 6 & 7
%
% For the single bits (0, 1, 2, and 3):
%
%     0 = No, this condition does not exist
%     1 = Yes, this condition exists.
%
% The double bits (6-7) represent aerosol levels w/ 4 values

% 0 - 00 Climatology
% 1 - 01 Low
% 2 - 10 Medium
% 3 - 11 High

flagname = {'fill','validAersol','water','interpolatedAerosol',...
    'aerosolLevel'};
datatype = {'logical','logical','logical','logical','uint8'};
nbits = [1 1 1 1 2];
bitPosition = [0 1 2 5 6];

N = LS8_AERO; % original 8-bit integers

I = 2.^nbits-1;
for k=1:length(flagname)
        %move bits to the right and only consider rightmost bits
        bitsToConsider = bitshift(N,-bitPosition(k));
        if strcmp(datatype{k},'logical')
            S.(flagname{k}) = bitand(bitsToConsider,I(k))>0;
        else
            S.(flagname{k}) = cast(bitand(bitsToConsider,I(k)),datatype{k});
        end
end

end