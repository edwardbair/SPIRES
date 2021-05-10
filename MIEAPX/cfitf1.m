function F1 = cfitf1( T )
%      Fits the function
%            F1( T ) = -exp(i*pi/6 ) * Ai-prime(Z1) / Ai(Z1)
%       where Ai is the Airy function,  Ai-prime its derivative
%           Z1 = T^2 * exp( 2*i*pi/3 )

Z1 = T.^2*exp(complex(0,2)*pi/3);
F1 = -exp(complex(0,1)*pi/6).*airy(1,Z1)./airy(Z1);

end