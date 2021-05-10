function [D] = D_calc(m,x,N)
%Compute the logarithmic derivative of the Ricatti-Bessel function at all orders (from 0 to N)

n = real(m);
kappa = abs(imag(m));
if n<1 || n>10 || kappa>10 || x*kappa>=3.9-10.8*n+13.78*n^2
    D = Ddownward(m*x,N);
else
    D = Dupward(m*x,N);
end
end

function D = Ddownward(z,N)
%returns the Ricatti-Bessel function of orders 0 to N for argument z using
%downwards recurrence relations

D = complex(zeros(N,1),zeros(N,1));
D(end) = LentzDn(z,N);
for n=N:-1:2
    D(n-1) = n/z-1/(D(n)+n/z);
end
end

function D = Dupward(z,N)
%returns the Ricatti-Bessel function of orders 0 to N for argument z using
%upwards recurrence relations

D = complex(zeros(N,1),zeros(N,1));
Dexp = exp(complex(0,-2)*z);
D(1) = -1/z+(1-Dexp)/((1-Dexp)/z-1i*(1+Dexp));
for n=2:N
    D(n) = 1/(n/z-D(n-1))-n/z;
end
end