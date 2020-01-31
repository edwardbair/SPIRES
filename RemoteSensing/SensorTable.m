function [ Tbl ] = SensorTable( sensor,varargin )
% [ Tbl ] = SensorTable( sensor,varargin )
%table of wavelength ranges and band designations for specified sensor
%
%Input
%   sensor - specify 'help' for list of available sensors
%Optional input
%   units - wavelength units, 'um' is default
%
%Output
%   Table showing band designations, wavelength ranges, and spatial
%   resolutions of specified sensor

narginchk(1,2)
nargoutchk(0,1)

p = inputParser;
addRequired(p,'sensor',@(x) ischar(x) || iscategorical(x));
defaultUnits = 'mum';
addOptional(p,'units',defaultUnits,@ischar)
parse(p,sensor,varargin{:});

if (ischar(p.Results.sensor) &&...
        contains(p.Results.sensor,'aviris','IgnoreCase',true)) ||...
        (iscategorical(p.Results.sensor) &&...
        contains(char(p.Results.sensor),'aviris','IgnoreCase',true))
    m = matfile('aviris-ng.mat');
    w = m.w;
    fwhm = m.fwhm;
    Band = (1:length(w))';
    Sensor = repmat({'AVIRIS-NG'},size(Band));
    Sensor = categorical(Sensor);
    CentralWavelength = w/1000;
    LowerWavelength = CentralWavelength-fwhm/(2*1000);
    UpperWavelength = CentralWavelength+fwhm/(2*1000);
    Bandwidth = fwhm/1000;
    SpatialResolution = repmat({'varies'},size(Band));
    SpatialResolution = categorical(SpatialResolution);
    Tbl = table(Sensor,Band,LowerWavelength,UpperWavelength,SpatialResolution,...
        CentralWavelength,Bandwidth);
    if ~strcmp(p.Results.units,defaultUnits)
        Tbl.CentralWavelength = convertUnits(Tbl.CentralWavelength,defaultUnits,p.Results.units);
        Tbl.LowerWavelength = convertUnits(Tbl.LowerWavelength,defaultUnits,p.Results.units);
        Tbl.UpperWavelength = convertUnits(Tbl.UpperWavelength,defaultUnits,p.Results.units);
        Tbl.Bandwidth = convertUnits(Tbl.Bandwidth,defaultUnits,p.Results.units);
        Tbl.Properties.VariableUnits = {'','',p.Results.units,p.Results.units,...
            'm',p.Results.units,p.Results.units};
    end
else
    m = matfile('multispectral.mat');
    MultispectralTable = m.MultispectralTable;
    availableSensors = unique(MultispectralTable.Sensor);
    if strcmpi(p.Results.sensor,'help')
        str = sprintf('\n%s\n',availableSensors(1));
        for k=2:length(availableSensors)
            str = cat(2,str,sprintf('%s\n',availableSensors(k)));
        end
        str = cat(2,str,sprintf('%s\n','AVIRIS-NG'));
        fprintf('%s: available sensors are: %s',...
            mfilename, str)
        return
    end
    
    if ischar(p.Results.sensor)
        t = strcmpi(cellstr(MultispectralTable.Sensor),p.Results.sensor);
    elseif iscategorical(p.Results.sensor)
        t = strcmpi(cellstr(MultispectralTable.Sensor),char(p.Results.sensor));
    end
    if nnz(t)==0
        if contains(char(p.Results.sensor),'fengyun','IgnoreCase',true)
            t = strcmpi(cellstr(MultispectralTable.Sensor),'fengyun-4');
        else
            error('sensor ''%s'' not recognized, call with ''help''',char(p.Results.sensor))
        end
    end
    Tbl = MultispectralTable(t,:);
    if ~strcmp(p.Results.units,defaultUnits)
        Tbl.LowerWavelength = convertUnits(Tbl.LowerWavelength,defaultUnits,p.Results.units);
        Tbl.UpperWavelength = convertUnits(Tbl.UpperWavelength,defaultUnits,p.Results.units);
        Tbl.Properties.VariableUnits = {'','',p.Results.units,p.Results.units,'m'};
    end
    
    % calculate the central wavelengths and bandwidth
    Tbl.CentralWavelength = mean([Tbl.LowerWavelength Tbl.UpperWavelength],2);
    Tbl.Bandwidth = Tbl.UpperWavelength-Tbl.LowerWavelength;
    Tbl.Properties.VariableUnits = {'','',p.Results.units,p.Results.units,'m',...
        p.Results.units,p.Results.units};
end
end