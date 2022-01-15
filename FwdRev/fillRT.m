function [newST] = fillRT(Results)
%fill the radiative transfer structure
newST.icer = Results.icer;
assert(~(istable(Results.lap) && Results.lookup),...
    'if ''LAP'' is specified as a table, ''lookup'' must be false')
newST.lookupMie = Results.lookup;

% calculation method
method = Results.method;
if strncmpi(method,'two',3) || strncmpi(method,'mea',3) || strncmpi(method,'hyb',3)
    newST.method = 'Meador-Weaver hybrid twostream';
elseif strncmpi(method,'del',3) || strncmpi(method,'edd',3)
    newST.method = 'delta-Eddington twostream';
elseif strncmpi(method,'dis',3)
    newST.method = 'disort';
else
    error('''method'' = ''%s'' unrecognized',method)
end

% diffuse fraction?
if isempty(Results.diffusefraction)
    newST.calcTerrain = false;
else
    newST.calcTerrain = true;
    newST.diffF = Results.diffusefraction;
end
end