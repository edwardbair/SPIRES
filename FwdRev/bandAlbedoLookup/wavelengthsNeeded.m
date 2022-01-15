function [wavelength] = wavelengthsNeeded(waveVector,bandPass)
%set just the wavelengths we need to cover the bands in bandPass
%Input
% waveVector - vector of wavelengths available
% bandPass - matrix of size(nBands,2) that specifies the bands (note that
%   they won't necessarily be in spectral order)
% (assumption: waveVector and bandPass are in the same units)
%Output
% wavelength - values from the original waveVector that cover the band
%   passes plus a little excess

p = inputParser;
addRequired(p,'waveVector',@(x) isnumeric(x) && isvector(x))
addRequired(p,'bandPass',@(x) isnumeric(x) && ismatrix(x) && size(x,2)==2)
parse(p,waveVector,bandPass)

%make sure waveVector covers the band passes
assert(min(waveVector)<=min(bandPass(:)),...
    'min(waveVector) = %g, is > min(bandPass) = %g (check units maybe)',...
    min(waveVector),min(bandPass(:)))
assert(max(waveVector)>=max(bandPass(:)),...
    'max(waveVector) = %g, is < max(bandPass) = %g (check units maybe)',...
    max(waveVector),max(bandPass(:)))

%cycle through the bands getting just what is needed
for b=1:size(bandPass,1)
    k1 = find(waveVector<min(bandPass(b,:)),1,'last');
    k2 = find(waveVector>max(bandPass(b,:)),1,'first');
    if isempty(k1)
        k1 = 1;
    elseif k1>1
        k1 = k1-1; % little extra at lower end of band
    end
    if isempty(k2)
        k2 = length(waveVector);
    elseif k2<length(waveVector)
        k2 = k2+1; % little extra and upper end of band
    end
    if b==1
        wavelength = waveVector(k1:k2);
    elseif isrow(waveVector)
        wavelength = cat(2,wavelength,waveVector(k1:k2));
    else
        wavelength = cat(1,wavelength,waveVector(k1:k2));
    end
end
wavelength = unique(wavelength); % elide duplicates and sort
end