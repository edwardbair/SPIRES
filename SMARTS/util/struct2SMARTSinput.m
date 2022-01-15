function nLines = struct2SMARTSinput(C,filename)
% nLines = struct2SMARTSinput(C,fileID)
%convert output from SetSMARTS into text file to send to smarts298.exe or
%smars295.exe
%C = output from SetSMARTS.m
%filename - full path to create input file for smarts298.exe or smarts295.exe
%nLines returns the number of lines written to the file

fileID = fopen(filename,'w');
fn = fieldnames(C);
N = length(fn);
nLines = 0;
for f=1:N
    mainV = C.(fn{f}).Value;
    isname = isfield(C.(fn{f}),'Name');
    if isname
        vname = C.(fn{f}).Name;
    end
    if ischar(mainV)
        fprintf(fileID,'''%s''  ',mainV);
    elseif isnumeric(mainV)
        fprintf(fileID,'%g ',mainV);
    elseif iscell(mainV)
        nVal = length(mainV);
        for k=1:nVal
            if isnumeric(mainV{k})
                fprintf(fileID,'%g ',mainV{k});
            elseif ischar(mainV{k})
                fprintf(fileID,'''%s'' ',mainV{k});
            end
        end
    end
    if isname
        fprintf(fileID,'\t\t! %s\n',vname);
    else
        fprintf(fileID,'\n');
    end
    nLines = nLines+1;
end
fclose(fileID);

% check for file EOF and fix if necessary
checkFileEOF(filename);
end