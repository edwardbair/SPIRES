function Tbl = getSMARTSspectrum(filename,varargin)
% Tbl = getSMARTSspectrum(filename,varargin)
%read spectrum from smarts298.ext.txt or smarts295.ext.txt file
%change column headers and output table
%if cosine(Z) is passed as the optional argument, TransDiffuse is also
%calculated if HorzDiffuse and TOA are available

% load the abbreviations
X = load('namesAbbrevUnits.mat');
abbrev = X.Tbl;

% read the SMARTS output spectrum, correcting some glitches
Tbl = readtable(filename,'PreserveVariablenames',true);
oldName = Tbl.Properties.VariableNames;
oldName = replace(oldName,'+','_');
oldName = replace(oldName,'-','_');
while any(contains(oldName,'__'))
    oldName = replace(oldName,'__','_');
end

% change the variable names
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
        if isundefined(abbrev{n,'Units'})
            units{k} = '';
        else
            units{k} = char(abbrev{n,'Units'});
        end
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