function [ rho ] = densityWater( TC )
% [ rho ] = densityWater( TC )
%density of water at atmospheric pressure, temp in deg C
%input can be scalar, vector, or matrix

% interpolated from table in Wikipedia
persistent already F

if isempty(already)
    X = [100 958.4; 80 971.8; 60 983.2; 40 992.2; 30 995.6502;...
        25 997.0479; 22 997.7735; 20 998.2071; 15 999.1026;...
        10 999.7026; 4 999.972; 0 999.8395; -10 998.117;...
        -20 993.547; -30 983.854];
    F = fit(X(:,1),X(:,2),'smoothingspline');
    already = true;
end

rho = F(TC);

end