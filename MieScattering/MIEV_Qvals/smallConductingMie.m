function [qext,qsca,qback,g] = smallConductingMie(x)
%Compute the efficiencies for small, perfectly conducting spheres

ahat1 = complex(0,(2/3)*(1-0.2*x^2))/complex(1-0.5*x^2,(2/3)*x^3);
bhat1 = complex(0,(x^2-10)/30)/complex(1+0.5*x^2,-x^3/3);
ahat2 = complex(0,x^2/30);
bhat2 = complex(0,-x^2/45);
qsca = x^4*(6*abs(ahat1)^2 + 6*abs(bhat1)^2+...
    10*abs(ahat2)^2 + 10*abs(bhat2)^2);
qext = qsca;
g = imag(ahat1)*(imag(ahat2)+imag(bhat1));
g = g + imag(bhat2)*((5/9)*imag(ahat2)+imag(bhat1));
g = g + real(ahat1)*real(bhat1);
g = g * (6*x^4/qsca);
qback = 9*x^4*abs(ahat1-bhat1-(5/3)*(ahat2-bhat2))^2;
end

