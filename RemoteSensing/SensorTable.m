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
        Tbl.CentralWavelength = convertLengthUnits(Tbl.CentralWavelength,...
            defaultUnits,p.Results.units);
        Tbl.LowerWavelength = convertLengthUnits(Tbl.LowerWavelength,...
            defaultUnits,p.Results.units);
        Tbl.UpperWavelength = convertLengthUnits(Tbl.UpperWavelength,...
            defaultUnits,p.Results.units);
        Tbl.Bandwidth = convertLengthUnits(Tbl.Bandwidth,defaultUnits,...
            p.Results.units);
    end
    Tbl.Properties.VariableUnits = {'','',p.Results.units,p.Results.units,...
        'm',p.Results.units,p.Results.units};
elseif (ischar(p.Results.sensor) &&...
        (contains(p.Results.sensor,'aso','IgnoreCase',true) ||...
        contains(p.Results.sensor,'casi','IgnoreCase',true))) ||...
        (iscategorical(p.Results.sensor) &&...
        contains(char(p.Results.sensor),'aso','IgnoreCase',true))
    m = matfile('aso-casi.mat');
    w = m.w;
    Band = (1:length(w))';
    Sensor = repmat({'ASO-CASI'},size(Band));
    Sensor = categorical(Sensor);
    CentralWavelength = w/1000;
    SpatialResolution = repmat({'varies'},size(Band));
    SpatialResolution = categorical(SpatialResolution);
    Tbl = table(Sensor,Band,SpatialResolution,...
        CentralWavelength);
    if ~strcmp(p.Results.units,defaultUnits)
        Tbl.CentralWavelength = convertLengthUnits(Tbl.CentralWavelength,...
            defaultUnits,p.Results.units);
    end
    Tbl.Properties.VariableUnits = {'','','m',p.Results.units};
elseif (ischar(p.Results.sensor) &&...
        contains(p.Results.sensor,'svc','IgnoreCase',true))
    m = matfile('svc-HR-1024i.mat');
    w = m.w;
    fwhm = m.fwhm;
    Band = (1:length(w))';
    Sensor = repmat({'SVC HR-1041i'},size(Band));
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
        Tbl.CentralWavelength = convertLengthUnits(Tbl.CentralWavelength,...
            defaultUnits,p.Results.units);
        Tbl.LowerWavelength = convertLengthUnits(Tbl.LowerWavelength,...
            defaultUnits,p.Results.units);
        Tbl.UpperWavelength = convertLengthUnits(Tbl.UpperWavelength,...
            defaultUnits,p.Results.units);
        Tbl.Bandwidth = convertLengthUnits(Tbl.Bandwidth,defaultUnits,...
            p.Results.units);
    end
    Tbl.Properties.VariableUnits = {'','',p.Results.units,p.Results.units,...
        'm',p.Results.units,p.Results.units};
elseif (ischar(p.Results.sensor) &&...
        (contains(p.Results.sensor,'enmap','IgnoreCase',true))) ||...
        (iscategorical(p.Results.sensor) &&...
        contains(char(p.Results.sensor),'enmap','IgnoreCase',true))
    vwidth = 6.5/1000;
    vmin = 420/1000;
    vmax = 1000/1000;
    lowV = (vmin:vwidth:vmax).';
    cV = lowV+vwidth/2;
    upV = lowV+vwidth;
    swidth = 10/1000;
    smin = 900/1000;
    smax = 2450/1000;
    lowS = (smin:swidth:smax).';
    cS = lowS+swidth;
    upS = lowS+swidth;
    for k=1:length(cV)
        x = num2str(k);
        vband(k,1) = categorical({[x 'VNIR']}); %#ok<AGROW>
    end
    for k=1:length(cS)
        x = num2str(k);
        sband(k,1) = categorical({[x 'SWIR']}); %#ok<AGROW>
    end
    SpatialResolution = repmat(30,length(vband)+length(sband),1);
    Sensor = categorical(repmat({'EnMAP'},length(vband)+length(sband),1));
    Band = cat(1,vband,sband);
    LowerWavelength = cat(1,lowV,lowS);
    UpperWavelength = cat(1,upV,upS);
    Bandwidth = cat(1,repmat(vwidth,length(vband),1),repmat(swidth,length(sband),1));
    CentralWavelength = cat(1,cV,cS);
    Tbl = table(Sensor,Band,LowerWavelength,UpperWavelength,SpatialResolution,...
        CentralWavelength,Bandwidth);
    if ~strcmp(p.Results.units,defaultUnits)
        Tbl.CentralWavelength = convertLengthUnits(Tbl.CentralWavelength,...
            defaultUnits,p.Results.units);
        Tbl.LowerWavelength = convertLengthUnits(Tbl.LowerWavelength,...
            defaultUnits,p.Results.units);
        Tbl.UpperWavelength = convertLengthUnits(Tbl.UpperWavelength,...
            defaultUnits,p.Results.units);
        Tbl.Bandwidth = convertLengthUnits(Tbl.Bandwidth,defaultUnits,...
            p.Results.units);
    end
else
    m = matfile('multispectral.mat');
    MultispectralTable = m.MultispectralTable;
    availableSensors = MultispectralTable.Sensor;
    availableSensors(end+1) = categorical({'AVIRIS-NG'});
    availableSensors(end+1) = categorical({'ASO-CASI'});
    availableSensors(end+1) = categorical({'EnMAP'});
    availableSensors = categorical(sort(unique(string(availableSensors))));
    if strcmpi(p.Results.sensor,'help')
        str = sprintf('\n%s\n',availableSensors(1));
        for k=2:length(availableSensors)
            str = cat(2,str,sprintf('%s\n',availableSensors(k)));
        end
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
            % next 2 lines need to be fixed after launch of Landsat 9
        elseif contains(char(p.Results.sensor),'landsat','IgnoreCase',true) &&...
                contains(char(p.Results.sensor),'oli','IgnoreCase',true)
            t = strcmpi(cellstr(MultispectralTable.Sensor),'landsat8oli');
        elseif contains(char(p.Results.sensor),'landsat','IgnoreCase',true) &&...
                contains(char(p.Results.sensor),'tirs','IgnoreCase',true)
            t = strcmpi(cellstr(MultispectralTable.Sensor),'landsat8tirs');
        else
            error('sensor ''%s'' not recognized, call with ''help''',...
                char(p.Results.sensor))
        end
    end
    Tbl = MultispectralTable(t,:);
    if ~strcmp(p.Results.units,defaultUnits)
        Tbl.LowerWavelength = convertLengthUnits(Tbl.LowerWavelength,...
            defaultUnits,p.Results.units);
        Tbl.UpperWavelength = convertLengthUnits(Tbl.UpperWavelength,...
            defaultUnits,p.Results.units);
        Tbl.Properties.VariableUnits = {'','',p.Results.units,p.Results.units,'m'};
    end
    
    % calculate the central wavelengths and bandwidth
    Tbl.CentralWavelength = mean([Tbl.LowerWavelength Tbl.UpperWavelength],2);
    Tbl.Bandwidth = Tbl.UpperWavelength-Tbl.LowerWavelength;
    Tbl.Properties.VariableUnits = {'','',p.Results.units,p.Results.units,'m',...
        p.Results.units,p.Results.units};
end
end