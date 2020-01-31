function [ g ] = getAsymmetry( x, qsca, an, bn )
% [ g ] = getAsymmetry( x, qsca, an, bn )
%calculates asymmetry parameter
%
% Bohren & Huffman, p 120, unnumbered equation

factr = 4/(x^2*qsca);
n = 1:length(an)-2;
anp = an;
anp(1) = [];
an(end-1:end) = [];
anp(end) = [];
bnp = bn;
bn(end-1:end) = [];
bnp(1) = [];
bnp(end) = [];
sum1v = sum((n.*(n+2)./(n+1)) .* real(an.*conj(anp)+ bn.*conj(bnp)));
sum2v = sum(((2*n+1)./(n.*(n+1))) .* real(an.*conj(bn)));

g = factr*(sum1v+sum2v);

end