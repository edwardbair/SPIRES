function [outStruct,varargout] = SMARTSMain(varargin)
% outStruct = SMARTSMain(varargin)
% outStruct = SMARTSMain(prescription,varargin)
% [outStruct,TblIntgRefl] = SMARTSMain(________)
%
%SMARTSMain runs SMARTS (Simple Model of Atmospheric Radiative Transfer of Sunshine)
%It works either for version 2.9.8 (preferred) or 2.9.5 (version names
%shortened to 298 or 295).
%For any session of MATLAB, call SetSMARTSversion first.
%
%SMARTS is a spectral solar radiation model developed by Christian Gueymard
%distributed by the DOE's National Renewable Energy Laboratory (NREL),
%https://www.nrel.gov/grid/solar-resource/smarts.html.
%It has a huge variety of input options, described in the Documentation
%folder of the folder where SNARTS (current version SMARTS_298_PC) is
%installed.
%To use this interface program, you must download SMARTS. It comes with
%an executable file smarts298.exe or smarts295.bat.exe.
%The Fortran source code for SMARTS is included in the download for 295, but
%not 298, so you
%could also compile it and call it anything you want. You need to put
%SMARTS in a folder that is writable by all users (sorry, a feature of
%SMARTS when running in batch mode).
%
%The documentation for the SMARTS input is briefly described in the
%accompanying function SetSMARTS.m (just type 'help SetSMARTS');
%the SMARTS user manual describes the input possibilities in detail. SMARTS
%also comes with an interactive .xls file that can be used to create an
%input file.
%
%Normally SMARTS runs from the input file smarts{298,295}.inp.ext, which contains
%the input variables in a prescribed order. This function takes the inputs
%as a set of Name/Value pairs in any order, creates the SMARTS input file,
%runs SMARTS via a system command, and retrieves the output as a structure.
%Included in the output structure is a prescription structure, which can
%also be used as the input to this function, and which can be modified for
%multiple runs of SMARTS with modest changes in the input.
%
%Inputs to this function
% Optionally, the next argument can be an input prescription structure,
%   perhaps modfied from the output of a previous run of SMARTS.
% Then the rest of the inputs consists of a set of Name/Value pairs to
% describe the SMARTS inputs, as documented in the SetSMARTS.m code and
% in the SMARTS user's manual.
%
%Output
% outStruct - structure than contains the output, including text and tables
%Optional output
% TblIntgRefl - input table for the SnowCloudReflectance toolbox for the
%   functions SnowCloudIntgRefl or invertSnowCloudIntgRefl
%
%IMPORTANT: The SMARTS input file and the output files from a run are in
%the smartsHome folder. When this function starts, it overwrites those
%files, so if you want to save them (you normally don't need to) you must
%rename them or move them to a different file. The files that are
%overwritten are: smarts{298,295}.inp.txt, smarts{298,295}.out.txt,
%smarts{298,295}.ext.txt, and smarts{298,295}.scn.txt.
%smarts298.inp.txt and smarts298.out.txt (or smarts295.inp.txt and smarts295.out.txt)
%are always produced, normally smarts{298,295}.ext.txt is produced, and
%smarts{298,295}.scn.txt might be depending on inputs.

%% which version
version = getSMARTSversion;
smartsHome = getSMARTShome;

%% The normal SMARTS files, delete if exist
switch version
    case 298
        inputSMARTS = 'smarts298.inp.txt';
        outputSMARTS = {'smarts298.out.txt','smarts298.ext.txt','smarts298.scn.txt'};
        scriptSMARTS = 'smartsexe.bat';
        programSMARTS = 'smarts298.exe';
    case 295
        inputSMARTS = 'smarts295.inp.txt';
        outputSMARTS = {'smarts295.out.txt','smarts295.ext.txt','smarts295.scn.txt'};
        scriptSMARTS = 'smartsexe.bat';
        programSMARTS = 'smarts295bat.exe';
end

% create the SMARTS script in a temporary file
[fileID,errmsg] = fopen(fullfile(smartsHome,scriptSMARTS),'w');
assert(fileID>0,errmsg)
if isunix
    fprintf(fileID,'cd "%s"\n./%s\n',smartsHome,programSMARTS);
else
    fprintf(fileID,'cd "%s"\n%s\n',smartsHome,programSMARTS);
end
smartsRun = fullfile(smartsHome,scriptSMARTS);
fclose(fileID);
if isunix
    system(['chmod a+x ' fullfile(smartsHome,scriptSMARTS)]);
end
%Input could include an input structure, and/or a set of Name/Value pairs
if isstruct(varargin{1})
    C = varargin{1};
    if length(varargin)>1
        % could fix this in future to allowing modification of prescription
        error('no additional arguments if you run with a prescription')
    end
else
    C = SetSMARTS(varargin{:});
end

% check if zenith angle is a field
iszang = false;
fn = fieldnames(C);
for k=1:length(fn)
    if isfield(C.(fn{k}),'Name')
        if contains(C.(fn{k}).Name,'ZENIT')
            iszang = true;
            zenith = C.(fn{k}).Value(1);
            break
        end
    end
end

%delete files as necessary
if exist(fullfile(smartsHome,inputSMARTS),'file')==2
    delete(fullfile(smartsHome,inputSMARTS));
end
for k=1:length(outputSMARTS)
    if exist(fullfile(smartsHome,outputSMARTS{k}),'file')==2
        delete(fullfile(smartsHome,outputSMARTS{k}));
    end
end

%% create the input file for SMARTS
nLines = struct2SMARTSinput(C,fullfile(smartsHome,inputSMARTS));
assert(nLines>0,'write to file ''%s'' failed',...
    fullfile(smartsHome,inputSMARTS))

%% run SMARTS
[status,cmdout] = system(smartsRun);
assert(status==0,cmdout)
delete(smartsRun);

%% retrieve output
% text output
if isunix
    tx = ['cat ' '"' fullfile(smartsHome,outputSMARTS{1}) '"'];
else
    tx = ['type ' '"' fullfile(smartsHome,outputSMARTS{1}) '"'];
end
[status,smartsText] = system(tx);
assert(status==0,smartsText);

% spectral table
if exist(fullfile(smartsHome,outputSMARTS{2}),'file')==2
    if iszang
        spectralTbl = getSMARTSspectrum(fullfile(smartsHome,...
            outputSMARTS{2}),cosd(zenith));
    else
        spectralTbl = getSMARTSspectrum(fullfile(smartsHome,outputSMARTS{2}));
    end
end

% scan table - still need to write getSMARTSscan function
if exist(fullfile(smartsHome,outputSMARTS{3}),'file')==2
    scanTbl = getSMARTSscan(fullfile(smartsHome,outputSMARTS{3}));
end

%% output structure
outStruct.prescription = C;
if exist('spectralTbl','var')==1
    outStruct.spectralTbl = spectralTbl;
end
if exist('scanTbl','var')==1
    outStruct.scanTbl = scanTbl;
end
outStruct.smartsText = smartsText;

% table for SnowCloudIntgRefl
if nargout>1
    assert(exist('spectralTbl','var')==1,'spectralTbl must exist for this output')
    Tbl = table(spectralTbl.waveL,[spectralTbl.HorzDirect spectralTbl.HorzDiffuse],...
        'VariableNames',{'wavelength','irradiance'});
    Tbl.Properties.VariableUnits = {'nm','W/m2/nm'};
    varargout{1} = Tbl;
end
end