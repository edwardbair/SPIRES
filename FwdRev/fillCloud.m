function [newSt] = fillCloud(st,Results,cloudType)
%fill cloud structure, all singles converted to doubles to support fmincon
assert(strcmpi(Results.sizeunit,'um') || strcmpi(Results.sizeunit,'mum'),...
    'unit for sizes of scatterers must be ''um'' or ''mum''')
newSt = st;
newSt.sizeUnit = Results.sizeunit;
newSt.radius = double(Results.radius);
newSt.temperature = double(Results.temperature);
if strcmpi(cloudType,'mixedCloud')
    newSt.wetness = double(Results.wetness);
    newSt.waterRadius = double(Results.waterradius);
end
assert(~isempty(Results.waterequivalent),...
    'for all cloud options, ''waterEquivalent'' must be specified')
newSt.waterEquivalent = double(Results.waterequivalent);
end