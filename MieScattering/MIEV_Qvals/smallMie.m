function [qext,qsca,qback,g] = smallMie(m,x)
%compute efficiencies for a small sphere with index m and size parameter x

m2 = m*m;
x2 = x*x;
D = m2+2+(1-0.7*m2)*x2-(8*m^4-385*m2+350)*x^4/1400+...
    complex(0,2)*(m2-1)*x^3*(1-0.1*x2)/3;
ahat1 = complex(0,2)*((m2-1)/3)*(1-0.1*x2+(4*m2+5)*x^4/1400)/D;
bhat1 = complex(0,1)*x2*((m2-1)/45)*(1+((2*m2-5)/70)*x2)/(1-(2*m2-5)/30*x2);
ahat2 = complex(0,1)*x2*((m2-1)/15)*(1-x2/14)/(2*m2+3-(2*m2-7)/14*x2);
T = abs(ahat1)^2+abs(bhat1)^2+5/3*abs(ahat2)^2;
temp = ahat2+bhat1;
g = real(ahat1*conj(temp))/T;
qsca = 6*x^4*T;
if imag(m)==0
    qext = qsca;
else
    qext = 6*x*real(ahat1+bhat1+5*ahat2/3);
end
sback = 1.5*x^3*(ahat1-bhat1-5*ahat2/3);
qback = 4*abs(sback)^2/x2;

end