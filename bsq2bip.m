function outImg = bsq2bip(inImg)
% outImg = bsq2bip(inImg)
%converts band-sequential image to band-interleaved-by-pixel
%
assert(ndims(inImg)==3,'input image must have 3 dimensions');
inSize = size(inImg);
if inSize(3)==1
    outImg = squeeze(inImg);
else
    outImg = permute(inImg,[3 1 2]);
end
end