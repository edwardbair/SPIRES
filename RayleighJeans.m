function f = RayleighJeans( lambda,T )
% f = RayleighJeans( lambda,T )
%Rayleigh-Jeans approximation to the Planck equation
%
% Wavelengths are in meters, temperatures in Kelvin.
% If the inputs are both scalars, result is a scalar.
% If one input is a vector or matrix and the other is a scalar,
% the result is a vector or matrix of the same size as the input.
% If both inputs are vectors or matrices, they must be the same size.
% If you want all combinations of 2 vectors, use meshgrid first.

% fundamental constants:
k=1.3806504e-23; % Boltzmann, J/K
c=299792458; % speed of light, m/s

% if both vectors or matrices, check sizes
[lambda,T] = checkSizes(lambda,T);
f = 2*c*k.*T./lambda.^4;

% error check
if nnz(lambda<=0) || nnz(isinf(lambda)) || nnz(T<=0)
    % check for zero or infinite wavelength, zero temperature, or negative values
    tz = lambda==0 | isinf(lambda) | T==0; % zero or infinite wavelength or zero temperature
    tn = lambda<0 | T<0; % negative wavelength or temperature
    
    % set result for zero or infinite wavelength or zero temperature to zero
    if nnz(tz)
        f(tz)=0;
    end
    % set results for negative wavelength or temperature to NaN
    if nnz(tn)
        warning('RayleighJeans: some wavelengths and/or temperatures are negative')
        f(tn)=NaN;
    end
end

end