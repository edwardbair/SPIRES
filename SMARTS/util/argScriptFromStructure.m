function argc = argScriptFromStructure(P)
%re-creates SMARTS298 or SMARTS295 input argument cell vector from structure
% takes advantage of the fact that the order doesn't matter

fn = fieldnames(P);
n = 0;
for k=1:length(fn)
    value = P.(fn{k}).Value;
    name = P.(fn{k}).Name;
    nV = length(value);
    if nV==1
        n = n+1;
        argc{n} = name; %#ok<AGROW>
        n = n+1;
            argc{n} = value; %#ok<AGROW>
    else
        %need to separate the variable names, but variable could be a
        %vector
        fb = strfind(strtrim(name),' ');
        if isempty(fb)
            n = n+1;
            argc{n} = name; %#ok<AGROW>
            n = n+1;
            argc{n} = value; %#ok<AGROW>
        else
            fb = cat(2,1,fb,length(name)); %start from left end
            for kV=1:nV
                n = n+1;
                argc{n} = strtrim(name(fb(kV):fb(kV+1)));
                n = n+1;
                argc(n) = {value(kV)};
            end
        end
    end
end
end