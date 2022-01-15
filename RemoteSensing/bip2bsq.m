function outImg = bip2bsq(inImg)
% outImg = bsq2bip(inImg)
%converts band-interleaved-by-pixel image to band-sequential
%
assert(ndims(inImg)==3,'input image must have 3 dimensions');
inSize = size(inImg);
if any(inSize(:)==1)
    outImg = squeeze(inImg);
    warning('image squeezed from size [%d %d %d] to [%d %s]',inSize,size(outImg)) 
else
    outImg = permute(inImg,[2 3 1]);
end
end