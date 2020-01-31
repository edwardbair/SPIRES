function [ S ] = unpackMOD09QC( QC_500m )
%unpack the bit fields in the MOD09
%   from MOD09 User's Guide v1.3, Table 15
%   http://modis-sr.ltdri.org/products/MOD09_UserGuide_v1_3.pdf
%
% Input
%   QC_500m - QC variable read from GetMOD09GA
%
% Output
%   S - structure with the following fields
%       modland QA - 0 good, 1 not so good, 2 or 3 bad
%       band QA (vector for 7 bands) - 0 best, 7 noisy detector, 8 dead
%           detector, 9 solar Z >= 86, 10 solar Z>=85 && <86, 11 missing
%           input, 12 climatology for atmosphere, 13 out of bounds, 14 data
%           faulty, 15 not processed, deep ocean or cloud
%       atmospheric correction - true or false
%       adjacency correction - true or false

flagname = {'modland','bandQA','atmosCorr','adjCorr'};
datatype = {'uint8','uint8','logical','logical'};
nbits = [2 4 1 1];

nbands = 7;
N = QC_500m; % original 32-bit integers
for k=1:length(flagname)
    I = 2^nbits(k)-1; % fill the lower nbits(k), leave others 0
    if k~=2
        if strcmp(datatype{k},'logical')
            S.(flagname{k}) = false(size(QC_500m));
            S.(flagname{k}) = bitand(N,I)>0;
        else
            S.(flagname{k}) = zeros(size(QC_500m),datatype{k});
            S.(flagname{k}) = cast(bitand(N,I),datatype{k});
        end
        N = bitshift(N,-nbits(k));
    else
        bq = zeros([size(QC_500m) nbands],datatype{k});
        for b=1:nbands
            bq(:,:,b) = cast(bitand(N,I),datatype{k});
            N = bitshift(N,-nbits(k));
        end
        S.(flagname{k}) = bq;
    end
end

end