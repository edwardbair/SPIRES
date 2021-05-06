function [An,Bn] = mieAnBn(m,x)
%compute Mie A,B coefficients at orders 0 to N
%
%Input
% m - complex index of refraction, imaginary part negative
% x - Mie size parameter (2*pi*r/lambda)
%
%output
% An,Bn - Mie A,B coefficients

nstop = floor(x+4.05*x^(1/3)+2)+1;
if real(m)>0
    D = D_calc(m,x,nstop+1);
end
An = complex(zeros(nstop,1),zeros(nstop,1));
Bn = complex(zeros(nstop,1),zeros(nstop,1));
psi0 = sin(x);
psi1 = psi0/x-cos(x);
xi0 = complex(psi0,cos(x));
xi1 = complex(psi1,cos(x)/x+sin(x));
for n=1:nstop
    if real(m)==0
        An(n) = (n*psi1/x-psi0)/(n*xi1/x-xi0);
        Bn(n) = psi1/xi1;
    elseif imag(m)==0
        z1r = real(D(n))/real(m)+n/x;
        An(n) = (z1r*psi1-psi0)/(z1r*xi1-xi0);
        z1r = real(D(n))*real(m)+n/x;
        Bn(n) = (z1r*psi1-psi0)/(z1r*xi1-xi0);
    else
        z1c = D(n)/m+n/x;
        An(n) = (z1c*psi1-psi0)/(z1c*xi1-xi0);
        z1c = D(n)*m+n/x;
        Bn(n) = (z1c*psi1-psi0)/(z1c*xi1-xi0);
    end
    xi = (2*n+1)*xi1/x-xi0;
    xi0 = xi1;
    xi1 = xi;
    psi0 = psi1;
    psi1 = real(xi1);
end
end