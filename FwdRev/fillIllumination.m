function [newST] = fillIllumination(fscript,Results)
%fill illumination geometry structure, all singles converted to doubles to
%support fmincon
newST.cosZ = double(Results.cosZ);
if fscript.substance==categorical({'snow'})
    if isempty(Results.muS)
        newST.muS = newST.cosZ;
    else
        newST.muS = double(Results.muS);
    end
    newST.viewF = double(Results.viewF);
else
    newST.muS = newST.cosZ;
    newST.viewF = 1;
end
newST.elevation = double(Results.elevation);
end