function [ S ] = unpackMOD09state( state_1km )
%unpack the bit fields in the MOD09 state_1km variable
%   from MOD09 User's Guide v1.3, Table 16
%   http://modis-sr.ltdri.org/products/MOD09_UserGuide_v1_3.pdf
%
% Input
%   state_1km - state variable read from GetMOD09GA
%
% Output
%   S - structure with the following fields
%       cloud - 0 clear, 1 cloudy, 2 mixed, 3 not sure
%       cloud shadow - true or false
%       landwater - 0 shallow ocean, 1 land, 2 coastline, 3 shallow lake,
%           4 ephemeral water, 5 deep lake, 6 continental ocean, 7 ocean
%       aerosol - 0 climatology, 1 low, 2 average, 3 high
%       cirrus - 0 none, 1 small, 2 average, 3 high
%       snow - from MOD35, true or false
%       (other bits not used)

flagname = {'cloud','cloudshadow','landwater','aerosol','cirrus','snow'};
datatype = {'uint8','logical','uint8','uint8','uint8','logical'};

nbits = [2 1 3 2 2 1];

for k=1:length(flagname)
    if strcmp(datatype{k},'logical')
        S.(flagname{k}) = false(size(state_1km));
    else
        S.(flagname{k}) = zeros(size(state_1km),datatype{k});
    end
end

N = state_1km; % original 16-bit integers
for k=1:length(flagname)
    I = 2^nbits(k)-1; % fill the lower nbits(k), leave others 0
    if strcmp(datatype{k},'logical')
        S.(flagname{k}) = false(size(state_1km));
        if strcmp(flagname{k},'cloudshadow')
            S.(flagname{k}) = bitand(N,I)>0;
        elseif strcmp(flagname{k},'snow') % MOD35 snow flag
            S.(flagname{k}) = bitget(state_1km,13)>0;
        else
            error('flagname %s not recognized as logical',flagname{k})
        end
    else
        S.(flagname{k}) = zeros(size(state_1km),datatype{k});
        S.(flagname{k}) = cast(bitand(N,I),datatype{k});
    end
    N = bitshift(N,-nbits(k));
end

end