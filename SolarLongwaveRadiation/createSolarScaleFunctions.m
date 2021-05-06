function [outS] = createSolarScaleFunctions(file,atmos)
% [outS] = createSolarScaleFunctions(file,atmos)
%create interpolation functions for SolarScale
%   Using the output from prepSolarScale, create the interpolating
%   functions
[~,~,ext] = fileparts(file);
assert(strcmpi(ext,'.mat'),'input file must be .mat')
m = matfile(file); % load the data to create the interpolating functions
inS = m.S;
FG = griddedInterpolant(inS.wavelength,inS.cosZ,inS.elev,inS.Global,'makima','nearest');
FD = griddedInterpolant(inS.wavelength,inS.cosZ,inS.elev,inS.DiffuseFraction,'makima','nearest');
Ttop = m.Ttop;
FTOA = griddedInterpolant(Ttop.wavelength,Ttop.TOA,'makima','nearest');
outS.Comment = ['functions for solar scale for Global and DiffuseFraction for '...
    upper(atmos) ' atmosphere and top-of-atmosphere'];
outS.FG = FG;
outS.FD = FD;
outS.FTOA = FTOA;
outS.waveUnits = 'nm';
outS.elevUnits = 'km';
end