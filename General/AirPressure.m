function P=AirPressure(P0,T0,z0,T,z)
% pressure at elevation
% P=AirPressure(P0,T0,z0,T,z)
%
% input
%   P0 (any units), T0 pressure and temperature (Kelvin) at elevation z0 (m)
%   T temperature(s) (Kelvin) at elevation z(s) (m)
%   if T is empty, use std lapse rate
%   z elevation at which pressure is needed, can be a vector
%   (if T also a vector, must be same size as z)

% values within range
if nnz(T0<=0) ~= 0 || (~isempty(T) && nnz(T<0) ~= 0)
    error('temperatures must be in degrees Kelvin')
end
%
% some constants
R=8.3143; % universal gas constant
M=0.02894; % average molecular weight of air
g=9.80665; % gravitational acceleration
%
%
% degenerate case, same elevation
if isscalar(z) && z==z0
    P=P0;
else
    %
    % ambient atmospheric lapse rate
    if isempty(T)
        gamma=-.0065;
        T = T0;
    elseif isscalar(T) % constant lapse rate
        gamma=(T-T0)/(z(1)-z0);
    else
        assert(isequal(size(T),size(z)),...
            'when T & z are both vectors, must be same size')
    end
    
    if isscalar(T)
        if gamma==0 % constant temperature
            T=T0;
            ratio = exp(-g*M*(z-z0)/(R*T));
        else % constant lapse rate
            ratio = (1+gamma*(z-z0)/T0).^(-g*M/(gamma*R));
        end
        P=P0.*ratio;
    else
        if isscalar(P0)
            P0 = repmat(P0,size(z));
        end
        if isscalar(T0)
            T0 = repmat(T0,size(T));
        end
        if isscalar(z0)
            z0 = repmat(z0,size(z));
        end
        ratio = zeros(size(z));
        for k=1:length(z)
            if (T(k)==T0(k))
                ratio(k) = exp(-g*M*(z(k)-z0(k))/(R*T(k)));
            else
                gamma = (T(k)-T0(k))/(z(k)-z0(k));
                ratio(k) = (1+gamma*(z(k)-z0(k))/T(k))^(-g*M/(gamma*R));
            end
        end
        P = P0.*ratio;
    end
end