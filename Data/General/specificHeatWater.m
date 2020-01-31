function [ Cwater ] = specificHeatWater( TC )
% [ Cwater ] = specificHeatWater( TC )
%specific heat of water, temp in deg C
%input can be scalar, vector, or matrix

% interpolated from table in Wikipedia which is in J/g/K
persistent already F

if isempty(already)
    X = [0 4.2176; 10 4.1921; 20 4.1818; 25 4.1814; 30 4.1784; 40 4.1785;...
        50 4.1806; 60 4.1843; 70 4.1895; 80 4.1963; 90 4.205; 100 4.2159];
    X(:,2) = X(:,2)*1000;
    F = fit(X(:,1),X(:,2),'smoothingspline');
    already = true;
end

Cwater = F(TC);

end