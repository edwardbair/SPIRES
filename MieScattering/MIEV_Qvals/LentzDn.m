function [Dn] = LentzDn(z,N)
%Compute the logarithmic derivative of the Ricatti-Bessel function of order
%N with argument z using the continued fraction technique of Lentz,
%Appl. Opt., 15, 668-671, (1976)

zinv = double(2)/z;
DN = double(N);
alpha = (DN+1/2)*zinv;
aj = -(DN+1+1/2)*zinv;
alpha_j1 = aj+1/alpha;
alpha_j2 = aj;
ratio = alpha_j1/alpha_j2;
runratio = alpha*ratio;

while abs(abs(ratio)-1)>1e-12
    aj = zinv-aj;
    alpha_j1 = 1/alpha_j1 + aj;
    alpha_j2 = 1/alpha_j2 + aj;
    ratio = alpha_j1/alpha_j2;
    zinv = -zinv;
    runratio = ratio*runratio;
end
Dn = double(-N)/z + runratio;
end