function outImg = bil2bsq(inImg)
% outImg = bsq2bip(inImg)
%converts band-interleaved-by-line image to band-sequential
%
assert(ndims(inImg)==3,'input image must have 3 dimensions');
inSize = size(inImg);
if inSize(2)==1
    outImg = squeeze(inImg);
else
    outImg = permute(inImg,[1 3 2]);
end
end