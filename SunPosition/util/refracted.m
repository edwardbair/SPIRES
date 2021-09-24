function muR = refracted( mu0, P, T )
% muR = refracted( mu0, P, T )
%   Adjust cosine of solar zenith angle for atmospheric refraction
%   Equation (42) from: Reda, I., and A. Andreas (2008), Solar position
%   algorithm for solar radiation applications (revised),
%   NREL/TP-560-34302, 56 pp, National Renewable Energy Laboratory, Golden,
%   CO, doi: 10.2172/15003974.
%   (in their equation, T is in deg C)
%
% Input:    mu0 - cosine of solar zenith angle, unrefracted
%           P - pressure in hPa (same as mb)
%           T - temperature in Kelvin
% (OK to use estimates for P & T as the sensitivity is small)
%
% Output: cosine of refracted solar zenith angle

% inputs can be scalars or vectors but, if more than one vector, they must
% be same size
t = isscalar(mu0)+isscalar(P)+isscalar(T);
if t==0
    assert(isequal(size(mu0),size(P),size(T)),...
        'all inputs vectors, so they must be the same size')
elseif t==1 % 2 of the 3 inputs are not scalar
    if isscalar(T)
        assert(isequal(size(mu0),size(P)),...
            'mu0 and P are both vectors but of different sizes')
    elseif isscalar(P)
        assert(isequal(size(mu0),size(T)),...
            'mu0 and T are both vectors but of different sizes')
    elseif isscalar(mu0)
        assert(isequal(size(T),size(P)),...
            'P and T are both vectors but of different sizes')
    end
end

% check P & T within range
tp = P>1200 | P<200;
assert(nnz(tp)==0, 'pressure must be in hPa, current range is [%g %g]',...
    min(P(:)), max(P(:)))
tt = T<30;
assert(nnz(tt)==0, 'Temp in Kelvin, current minimum value is %g', min(T(:)))

% Code following takes advantage of the sin/cos relationship,
% i.e. sind(90-x) = cosd(x) etc.

% solar elevation, degrees
e0 = asind(mu0);

% refraction - equation (42) from Reda-Andreas
delE = (P/1010).*(283./T).*(1.02./(60*tand(e0 + 10.3./(e0+5.11))));

% cosine of refracted solar zenith angle
muR = sind(e0+delE);

end