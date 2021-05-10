function [ value ] = convertUnits( inputValue, unitsFrom, varargin )
% [ value ] = convertUnits( inputValue, unitsFrom [, unitsTo] )
%FUNCTION DEPRECATED because of possible conflict with existing MATLAB
%function, use convertLengthUnits instead
%
%convert input values to different units
%(based on MATLAB function unitsratio, but considers the range of metric
%units in radiation)
%
% Input
%   inputValue - scalar, vector, or matrix
%   unitsFrom and unitsTo - length unit, possible values are
%       angstrom
%       nm or nanometer
%       um or micron or micrometer or mum
%       mm or millimeter
%       cm or centimeter
%       km or kilometer
%       GHz or gigahertz
%       MHz or megahertz
%       THz or terahertz
%       invcm or cmm1
% Optional input
%   unitsTo (default m)
%
% Output
%   value - input Values converted to unitsTo (default m)

% uses the MATLAB function unitsratio, with the following modifications
%   add angstrom as a unit
%   nm = nanometer, not nautical mile
%   um or micrometer = micron
%   add GHz, MHz, and 
%   add invcm

warning('Function %s DEPRECATED because of possible conflict with existing MATLAB function, use convertLengthUnits instead',...
    mfilename);

if nargin>2
    unitsTo = varargin{1};
else
    unitsTo = 'm';
end

if strcmpi(unitsTo,unitsFrom)
    value = inputValue;
    return
end

c=299792458; % speed of light, m/s, needed for frequency conversions
switch unitsFrom
    case 'angstrom'
        unitsFrom = 'm';
        inputValue = inputValue*1.e-10;
    case {'nm','nanometer'}
        unitsFrom = 'm';
        inputValue = inputValue*1.e-9;
    case {'um','mum','micrometer'}
        unitsFrom = 'micron';
    case 'MHz'
        unitsFrom = 'm';
        inputValue = c./(inputValue*1.e6);
    case 'GHz'
        unitsFrom = 'm';
        inputValue = c./(inputValue*1.e9);
    case 'THz'
        unitsFrom = 'm';
        inputValue = c./(inputValue*1.e12);
    case {'invcm','cmm1'}
        unitsFrom = 'cm';
        inputValue = 1./inputValue;
end
frequencyConversion = false;
reciprocalConversion = false;
holdUnits = unitsTo;
conversion = 1;
switch unitsTo
    case 'angstrom'
        unitsTo = 'm';
        conversion = 1.e10;
    case {'nm','nanometer'}
        unitsTo = 'm';
        conversion = 1.e9;
    case {'um','mum','micrometer'}
        unitsTo = 'micron';
    case 'MHz'
        unitsTo = 'm';
        frequencyConversion = true;
    case 'GHz'
        unitsTo = 'm';
        frequencyConversion = true;
    case 'THz'
        unitsTo = 'm';
        frequencyConversion = true;
    case {'invcm','cmm1'}
        unitsTo = 'cm';
        reciprocalConversion = true;
end

value = inputValue*unitsratio(unitsTo,unitsFrom);
if frequencyConversion
    switch holdUnits
        case 'MHz'
            conversion = 1.e-6;
        case 'GHz'
            conversion = 1.e-9;
        case 'THz'
            conversion = 1.e-12;
    end
    value = conversion*c./value;
elseif reciprocalConversion
    value = 1./value;
else
    value = conversion*value;
end

end