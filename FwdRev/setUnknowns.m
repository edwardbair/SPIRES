function [solveFor,newPrescription] = setUnknowns(substance,unknowns,Prescription)
%validate the choice of unknowns for inversion
%
%Input
% substance - categorical: snow, iceCloud, waterCloud, or mixedCloud
% unknowns - cell vector of proposed unknowns
% Spectrum - structure from fwdPrescription, fields include wavelength and
%   waveUnit
%
%Output
% solveFor - scrubbed, validated list of unknowns
% newPrescription - perhaps modified from input

p = inputParser;
addRequired(p,'substance',@iscategorical)
addRequired(p,'unknowns',@iscell)
addRequired(p,'Prescription',@isstruct)
parse(p,substance,unknowns,Prescription)
Spectrum = Prescription.Spectrum;

% possible unknowns, varies with substance
switch substance
    case categorical({'snow'})
        validStrings = {'fSCA','wetness','LAPfraction',...
            'LAPradius','waterRadius','radius',...
            'waterEquivalent','fractionalCoverage','muS'};
        equivStrings = {'fSCA','fractionalCoverage'};
    case {categorical({'iceCloud'}),categorical({'waterCloud'})}
        validStrings = {'radius','waterEquivalent'};
    case categorical({'mixedCloud'})
        validStrings = {'radius','waterRadius','wetness','waterEquivalent'};
    otherwise
        error('substance ''%s'' not recognized, valid are ''snow'', ''iceCloud'', ''waterCloud'', or ''mixedCloud''',...
            char(substance))
end

% find valid unknowns
solveFor = cell(1,length(unknowns));
ns = 1; % count is adaptive because multiple laps will add extra variable
for k=1:length(unknowns)
    matchedstr = validatestring(unknowns{k},validStrings);
    if exist('equivStrings','var') &&...
            contains(matchedstr,equivStrings(:),'IgnoreCase',true)
        solveFor{ns} = equivStrings{1};
    elseif strcmpi(matchedstr,'lapfraction') || strcmpi(matchedstr,'lapradius')
        if Prescription.snow.nLAP == 1
            solveFor{ns} = matchedstr;
        elseif Prescription.snow.nLAP == 2
            solveFor{ns} = [matchedstr '1'];
            ns = ns+1;
            solveFor{ns} = [matchedstr '2'];
        else
            error('cannot reliably solve for more than 2 %s',matchedstr)
        end
    else
        solveFor{ns} = matchedstr;
    end
    ns = ns+1;
end
% if fSCA is a variable, it could turn into a vector, so put it at the end
t = contains(solveFor,'fSCA','IgnoreCase',true);
if any(t)
    solveFor = cat(2,solveFor(~t),{'fSCA'});
end

% make sure wavelenghts long enough to get size of scatterer, if asked for
wavelimitGrainSize = 1045; % nm
if any(contains(solveFor,'radius','IgnoreCase',true)) ||...
        any(contains(solveFor,'waterRadius','IgnoreCase',true))
    assert(convertLengthUnits(max(Spectrum.wavelength(:)),...
        Spectrum.waveUnit,'nm')>=wavelimitGrainSize,...
        'maximum wavelength must be >= %f %s to retrieve size of snow or cloud scatterer',...
        convertLengthUnits(wavelimitGrainSize,'nm',Spectrum.waveUnit),Spectrum.waveUnit)
end

% make sure wavelengths short enough to characterize LAP
wavelimitLAP = 500;
if any(contains(solveFor,'lap','IgnoreCase',true))
    assert(convertLengthUnits(min(Spectrum.wavelength(:)),...
        Spectrum.waveUnit,'nm')<=wavelimitLAP,...
        'minimum wavelength must be <= %f %s to retrieve LAP properties',...
        convertLengthUnits(wavelimitLAP,'nm',Spectrum.waveUnit),Spectrum.waveUnit)
end
newPrescription = Prescription;

% warn about solving for LAPfraction and LAPradius
if any(contains(solveFor,'lapfraction','IgnoreCase',true)) &&...
        any(contains(solveFor,'lapradius','IgnoreCase',true))
    warning('usually you cannot solve for ''LAPfraction'' AND ''LAPradius'' because different combinations can produce too similar results')
end
end