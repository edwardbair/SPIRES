function [Tbl,varargout] = LinearMixture_spectral(too,fractions,varargin)
%  Tbl = LinearMixture_spectral(TablesOrObjects,fractions [,'wavelength',value,'waveUnit',value])
%  [Tbl,F] = LinearMixture_spectral(TablesOrObjects,fractions [,'wavelength',value,'waveUnit',value])
%linear mixture of spectrometer data
%
%Calculate linear mixtures of spectral reflectances
%
%Input
% TablesOrObjects - specifying wavelengths and reflectances of the
%   endmembers, which can be:
%   a single table with 'wavelength' and 'reflectance', where reflectances
%       are specified as a matrix of nWavelength x nEndMember
%   a cell vector of tables, with the first table specifying the master set
%       of wavelengths
%   a cell vector of objects, either griddedInterpolants or cfits which
%       specify the reflectances as functions of wavelengths
% fractions - numerical vector of fractional values that sum to 1.0 (and
%   adjusted to that sum)
%Optional inputs, name/value pairs to specify either 'wavelength' or 'sensor'
%   just one can be specified
% 'wavelength' – vector
% 'waveUnit' – units for wavelength, no default (to prevent errors)
% 'sensor' – instead of 'wavelength', can specify a spectrometer
%    in the SensorTable.m function
%
%Output
% Tbl - table of wavelengths and reflectances of linear mixture
%Optional output
% F - griddedInterpolant object of reflectances of linear mixture as a
%   function of wavelength

p = inputParser;
addRequired(p,'too',@(x) istable(x) || iscell(x))
addRequired(p,'fractions',@isnumeric)
addParameter(p,validatestring('wavelength',{'wavel','wavelength'}),...
    [],@(x) isempty(x) || (isnumeric(x) && all(x>0)))
addParameter(p,validatestring('waveunit',{'waveu','waveunit','lambda'}),...
    '',@ischar)
addParameter(p,validatestring('sensor',{'sens','sensor'}),...
    '',@(x) ischar(x) || iscategorical(x))
parse(p,too,fractions,varargin{:})

% wavelengths from table or from optional argument
assert(isempty(p.Results.wavelength) || isempty(p.Results.sensor),...
    '''wavelength'' or ''sensor'' can be specified, but not both')
if xor(isempty(p.Results.wavelength),isempty(p.Results.sensor))
    assert(~isempty(p.Results.waveunit),...
        'if ''wavelength'' or ''sensor'' is specified, ''waveUnit'' must also be specified')
    waveUnit = p.Results.waveunit;
    if isempty(p.Results.wavelength)
        sTbl = SensorTable(p.Results.sensor,waveUnit);
        wavelength = sTbl.CentralWavelength;
    else
        wavelength = p.Results.wavelength;
    end
end
if istable(p.Results.too) % all reflectances in one table
    assert(size(too.reflectance,2)==length(fractions),...
        'if reflectances are specified in a single table, must have one column for each endmember fraction')
    if exist('waveUnit','var')
        Frefl = table2Interpolant(p.Results.too,waveUnit);
    else
        [Frefl,waveUnit,wavelength] = table2Interpolant(p.Results.too);
    end
else % cell vector of either tables or objects
    if istable(p.Results.too{1}) % cell vector of tables
        assert(length(p.Results.too)==length(fractions),...
            'if reflectances are specified in multiple tables, must have one table for each endmember fraction')
        % build gridded interpolant for reflectance of each endmember
        Frefl = cell(1,length(p.Results.too));
        for k=1:length(p.Results.too)
            thisTbl = p.Results.too{k};
            if exist('waveUnit','var')
                Frefl{k} = table2Interpolant(thisTbl,waveUnit);
            elseif k==1
                [Frefl{k},waveUnit,wavelength] = table2Interpolant(thisTbl);
            end
        end
    else % cell vector of objects, either griddedInterpolant or cfit
        assert(length(p.Results.too)==length(fractions),...
            'if reflectances are specified in multiple interpolating objects, must have one object for each endmember fraction')
        Frefl = p.Results.too;
        assert(~isempty(wavelength) && ~isempty(waveUnit),...
            'if input interpolant specified, so must ''wavelength'' and ''waveUnit''')
    end
end
% normalize fractions to sum to 1.0
sumF = sum(fractions);
fractions = fractions/sumF;
% calculate the results for each wavelength
if exist('wavelength','var') && ~isempty(wavelength)
reflectance = zeros(length(wavelength),1);
else
end
for f=1:length(fractions)
    thisF = Frefl{f};
    reflectance = reflectance+fractions(f)*thisF(wavelength);
end
Tbl = table(wavelength,reflectance);
Tbl.Properties.VariableUnits = {waveUnit,''};
if nargout>1 % interpolating object if desired
    F = griddedInterpolant(wavelength,reflectance,'makima','nearest');
    varargout{1} = F;
end
end

function [gI,varargout] = table2Interpolant(reflTbl,varargin)
p = inputParser;
addRequired(p,'reflTbl',@istable);
addOptional(p,'waveunit','',@ischar);
parse(p,reflTbl,varargin{:});
Frefl = cell(1,size(reflTbl.reflectance,2));
if isempty(p.Results.waveunit)
    for k=1:size(reflTbl.reflectance,2)
        Frefl{k} = griddedInterpolant(reflTbl.wavelength,reflTbl.reflectance(:,k),...
            'makima','nearest');
    end
    t = contains(reflTbl.Properties.VariableNames,'wavelength');
    if nargout>1
    varargout{1} = reflTbl.Properties.VariableUnits{t};
    end
    if nargout>2
        varargout{2} = reflTbl.wavelength;
    end
else
    t = contains(reflTbl.Properties.VariableNames,'wavelength');
    assert(nnz(t)==1,'one column of input table must be ''wavelength''')
    tblWave = convertLengthUnits(reflTbl.wavelength,reflTbl.Properties.VariableUnits{t},p.Results.waveunit);
    for k=1:size(reflTbl.reflectance,2)
        Frefl{k} = griddedInterpolant(tblWave,reflTbl.reflectance(:,k),...
            'makima','nearest');
    end
    if nargout>1
    varargout{1} = p.Results.waveunit;
    end
    if nargout>2
        varargout{2} = tblWave;
    end
end
% just one interpolant if table contains just one reflectance
if length(Frefl)==1
    gI = Frefl{1};
else
    gI = Frefl;
end
end