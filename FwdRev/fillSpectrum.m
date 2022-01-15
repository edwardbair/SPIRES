function [newST] = fillSpectrum(Results)
%fill wavelength / sensor structure, singles converted to doubles to
%accommodate fmincon
assert(xor(isempty(Results.wavelength),isempty(Results.sensor)) ||...
    xor(isempty(Results.wavelength),isempty(Results.bandpass)),...
    'must specify just one of ''wavelength'', ''sensor'', or ''bandPass''')

newST.waveUnit = Results.waveunit;
assert(~isempty(newST.waveUnit),'''waveUnit'' must be specified')
if ~isempty(Results.wavelength) && isvector(Results.wavelength)
    newST.spectrometer = true;
    newST.wavelength = double(Results.wavelength(:));
elseif ~isempty(Results.sensor)
    T = SensorTable(Results.sensor,newST.waveUnit);
    % other spectrometers to be added here, e.g. EnMAP, CHIME, SBG
    if contains(Results.sensor,'aviris','IgnoreCase',true) ||...
            contains(Results.sensor,'casi','IgnoreCase',true) ||...
            contains(Results.sensor,'svc','IgnoreCase',true)
        newST.spectrometer = true;
        newST.wavelength = double(T.CentralWavelength(:));
        newST.sensor = Results.sensor;
    else
        newST.spectrometer = false;
        newST.sensor = Results.sensor;
        if isempty(Results.bands)
            newST.bandPass = double([T.LowerWavelength T.UpperWavelength]);
            newST.bands = T.Band.';
        else
            Band = categorical(Results.bands);
            holdT = innerjoin(table(Band(:),'VariableNames',{'Band'}),T,'Keys',{'Band'});
            newST.bandPass = [holdT.LowerWavelength holdT.UpperWavelength];
            newST.bands = holdT.Band(:).';
        end
    end
elseif ~isempty(Results.bandpass)
    newST.spectrometer = false;
    newST.bandPass = double(Results.bandpass);
    newST.bands = 1:size(newST.bandPass,1);
end
% convert multispectral to pseudo-spectrometer by picking wavelengths
% across each band
if ~newST.spectrometer
    minWave = 3;
    delWaveThresh = 15; % nm
    for k=1:size(newST.bandPass,1)
        % set number of wavelengths for this band
        bandWidth = convertLengthUnits(newST.bandPass(k,2)-newST.bandPass(k,1),...
            newST.waveUnit,'nm');
        nWave = max(minWave,round(bandWidth/delWaveThresh));
        thisWave = linspace(newST.bandPass(k,1),newST.bandPass(k,2),nWave);
        if k==1
            wavelength = thisWave(:);
            wbands = repmat(newST.bands(k),1,nWave);
        else
            wavelength = cat(1,wavelength,thisWave(:));
            wbands = cat(2,wbands,repmat(newST.bands(k),1,nWave));
        end
    end
    newST.wavelength = wavelength;
    newST.waveBand = wbands.';
end
    
newST.ignoreSolar = logical(Results.ignoresolar);
end