function F2 = cfitf2( T )
%       Fits the function
%            F2( T ) =  - exp( i*pi/6 ) * Ai-prime(Z2) / Ai(Z2)
%      Ai is the Airy function,  Ai-prime its
%      derivative,  and
%            Z2 = T^2 * exp( - i*pi/3 )

Z2 = T.^2*exp(complex(0,-1)*pi/3);
F2 = -exp(complex(0,1)*pi/6).*airy(1,Z2)./airy(Z2);

end