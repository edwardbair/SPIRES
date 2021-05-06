function [FG,FD,FTOA,waveUnits,elevUnits] = getSolarScaleInterpolants(atmos)
%retrieve solar scale functions created by createSolarScaleFunctions
%atmosphere code is embedded in file name
%
file = ['SolarScale_' upper(atmos) '.mat'];
m = matfile(file);
FG = m.FG;
FD = m.FD;
FTOA = m.FTOA;
waveUnits = m.waveUnits;
elevUnits = m.elevUnits;
end

