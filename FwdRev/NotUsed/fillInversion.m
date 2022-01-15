function newST = fillInversion(prescription,Results)
%fill the inversion structure
Rfun = Results.Rfun;
Spectrum = prescription.Spectrum;
if Spectrum.Spectrometer
    newST.reflectance = Rfun(Spectrum.wavelength);
else
    newST = Rfun(Spectrum.bandPass);
end
% solution method
if contains(Results.solutionmethod,'norm')
    newST.solutionMethod = 'normResiduals';
elseif contains(Results.solutionmethod,'angle')
    newST.solutionMethod = 'spectralAngle';
elseif contains(Results.solutionmethod,'hyb')
    newST.solutionMethod = 'hybrid';
else
    error('solution method %s not recognized',Results.solutionmethod);
end
% possible unknowns, check compatible with prescription
if isfield(prescription,'snow')
    pu = {'radius','fSCA','wetness','waterEquivalent','dust',...
        'dustRadius','soot','sootRadius'};
elseif isfield(prescription,'iceCloud') || isfield(prescription,'waterCloud')
    pu = {'radius','waterEquivalent'};
elseif isfield(prescription,'mixedCloud')
    pu = {'radius','waterRadius','wetness','waterEquivalent'};
else
    error('no recognized field in prescription')
end
% unknowns could include 'nRwt' if hybrid solution method
if strcmp(newST.solutionMethod,'hybrid')
    pu = cat(2,pu,'nRwt');
end

% list the unknowns in the inversion output structure
assert(iscell(p.Results.unknowns),'''unknowns'' must be a cell vector')
for k=1:length(Results.unknowns)
    if k==1
        newST.unknowns = {validatestring(Results.unknowns{k},pu)};
    else
        newST.unknowns = cat(2,newST.unknowns,validatestring(Results.unknowns{k},pu));
    end
end

end