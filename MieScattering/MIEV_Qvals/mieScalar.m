function [qext,qsca,qback,g] = mieScalar(m,x)
%Mie terms for scalar argument
%translated from Scott Prahl's python version
%
%Input
% m - complex index of refraction, imaginary part negative
% x - Mie size parameter (2*pi*r/lambda)
%
%Output
% qext - Mie extinction efficiency
% qsca - Mie scattering efficiency
% qback - Mie backscattering efficiency
% g - Mie asymmetry parameter

if real(m)==0 && x<0.1
    [qext,qsca,qback,g] = smallConductingMie(x);
elseif real(m)>0 && abs(m)*x<0.1
    [qext,qsca,qback,g] = smallMie(m,x);
else
    [an,bn] = mieAnBn(m,x);
    nmax = length(an);
    n = (1:nmax).';
    cn = 2*n+1;
    qext = (2/x^2)*sum(cn.*(real(an)+real(bn)));
    if imag(m)==0
        qsca = qext;
    else
        qsca = (2/x^2)*sum(cn.*(abs(an).^2+abs(bn).^2));
    end
    qback = abs(sum((-1).^n.*cn.*(an-bn))).^2/x^2;
    c1n = n.*(n+2)./(n+1);
    c2n = cn./n./(n+1);
    g = 0;
    for k=1:nmax-1
        asy1 = c1n(k)*real(an(k)*conj(an(k+1))+bn(k)*conj(bn(k+1)));
        asy2 = c2n(k)*real(an(k)*conj(bn(k)));
        g = g+4*(asy1+asy2)/qsca/x^2;
    end
end
end