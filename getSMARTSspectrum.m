function Tbl = getSMARTSspectrum(filename,varargin)
% Tbl = getSMARTSspectrum(filename,varargin)
%read spectrum from smarts295.ext.txt file
%change column headers and output table
%if cosine(Z) is passed as the optional argument, TransDiffuse is also
%calculated if HorzDiffuse and TOA are available

persistent abbrev
% load the abbreviation table
if isempty(abbrev)
    X = load('namesAbbrevUnits.mat');
    abbrev = X.Tbl;
end

% read the SMARTS output spectrum
Tbl = readtable(filename);

% change the variable names
oldName = Tbl.Properties.VariableNames;
newName = cell(size(oldName));
units = cell(size(oldName));
for k=1:length(oldName)
    n = find(strcmp(abbrev.origName,oldName{k}));
    if isempty(n)
        newName{k} = oldName{k};
        warning('variable name %s not found in abbreviation table',oldName{k})
        units{k} = 'unknown';
    else
        newName{k} = char(abbrev{n,'tblName'});
        units{k} = char(abbrev{n,'Units'});
    end
end

Tbl.Properties.VariableNames = newName;
Tbl.Properties.VariableUnits = units;
Tbl.Properties.VariableDescriptions = oldName;

% diffuse transmittance?
if nargin>1
    cosZ = varargin{1};
    if cosZ<=0 || cosZ>1
        warning('cosZ = %f out of range',cosZ)
    elseif ~any(strcmpi(Tbl.Properties.VariableNames,'TOA')) ||...
            ~any(strcmpi(Tbl.Properties.VariableNames,'HorzDiffuse'))
        warning('flux variables to compute TransDiffuse not in table')
    else
        TransDiffuse = Tbl.HorzDiffuse./(cosZ*Tbl.TOA);
        addTable = table(TransDiffuse);
        addTable.Properties.VariableUnits = {'dimensionless'};
        addTable.Properties.VariableDescriptions = {'Diffuse_rad_transmittance'};
        Tbl = [Tbl addTable];
        Tbl.Properties.UserData = {'zenith' acosd(cosZ)};
    end
end

end