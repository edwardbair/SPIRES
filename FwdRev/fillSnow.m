function [newSt] = fillSnow(fscript,snowlim,Results)
%fill snow structure, all singles converted to doubles to support fmincon
assert(strcmpi(Results.sizeunit,'um') || strcmpi(Results.sizeunit,'mum'),...
    'unit for sizes of scatterers must be ''um'' or ''mum''')
newSt = fscript.snow;
newSt.sizeUnit = Results.sizeunit;
newSt.radius = double(Results.radius);
if isempty(Results.lap)
    newSt.cleanSnow = true;
    assert(isempty(Results.lapfraction),...
        '''LAP'' is empty but ''LAPfraction'' is not, must specify LAP')
else
    newSt.cleanSnow = false;
    if ischar(Results.lap)
        newSt.LAP = {Results.lap};
        newSt.nLAP = 1;
    elseif iscell(Results.lap)
        newSt.LAP = Results.lap;
        newSt.nLAP = length(newSt.LAP);
    elseif istable(Results.lap) % convert to matrix, rows along wavelength, columns for each LAP
        tmpT = Results.lap;
        assert(any(contains(tmpT.Properties.VariableNames,'wavelength')) &&...
            any(contains(tmpT.Properties.VariableNames,'RefractiveIndex')),...
            'table of refractive indices must have column labels ''wavelength'' and ''RefractiveIndex''')
        newSt.nLAP = size(tmpT.RefractiveIndex,2);
        cz = zeros(length(fscript.Spectrum.wavelength),newSt.nLAP);
        refidx = complex(cz,cz);
        tmpT = Results.lap;
        for k=1:newSt.nLAP
            Freal = griddedInterpolant(tmpT.wavelength,...
                real(tmpT.RefractiveIndex(:,k)),'pchip','nearest');
            Fimag = griddedInterpolant(tmpT.wavelength,...
                abs(imag(tmpT.RefractiveIndex(:,k))),'pchip','nearest');
            refidx(:,k) = complex(Freal(fscript.Spectrum.wavelength),...
                Fimag(fscript.Spectrum.wavelength));
        end
        newSt.LAP = refidx;
    end
    newSt.LAPfraction = double(Results.lapfraction);
    assert(newSt.nLAP==length(newSt.LAPfraction),...
        'number of LAPs must equal number of LAPfractions')
    if isempty(Results.lapradius)
        assert(~istable(newSt.LAP),...
            'if ''LAP'' is specified as a refractive index, ''LAPradius'' must be specified')
        for k=1:length(newSt.LAPfraction)
            switch newSt.LAP{k}
                case 'dust'
                    newSt.LAPradius(k) = snowlim.defaultDustRadius;
                case 'soot'
                    newSt.LAPradius(k) = snowlim.defaultSootRadius;
                otherwise
                    error('LAP %s not recognized',newSt.LAP(k))
            end
        end
    else
        newSt.LAPradius = double(Results.lapradius);
        assert(length(newSt.LAPradius)==newSt.nLAP,...
            'number of LAPs must equal length of LAPradius vector')
    end
end
newSt.temperature = Results.temperature;
if Results.wetness==0
    newSt.wetSnow = false;
else
    newSt.wetSnow = true;
    newSt.wetness = double(Results.wetness);
end
newSt.waterEquivalent = double(Results.waterequivalent);
if isinf(newSt.waterEquivalent)
    newSt.deepSnow = true;
else
    newSt.deepSnow = false;
end
% fractional snow
Substrate = fscript.Substrate;
if isempty(Substrate)
    if Results.fSCA(1)<1
        warning('''fSCA''<1, but Substrate substruct is empty, resetting to 1')
        newSt.fSCA = [1 0];
    end
    
    % substrate reflectance constants, check size of fSCA
else
    lenR0 = length(Substrate.fR0);
    assert(~(length(Results.fSCA)>1+lenR0),'length(''fSCA'') > 1 + length(''R0'')')
    if length(Results.fSCA)==1+lenR0
        newSt.fSCA = double(Results.fSCA(:).'); % row vector
    elseif length(Results.fSCA)==lenR0 % need 1 extra endmember
        tmpSCA = double(Results.fSCA(:).');
        % add endmember and sum to 1
        tmpSCA = [tmpSCA 1-sum(tmpSCA)];
        newSt.fSCA = tmpSCA/sum(tmpSCA);
    elseif length(Results.fSCA)<lenR0
        tmpSCA = zeros(1,1+lenR0);
        for k=1:length(Results.fSCA)
            tmpSCA(k) = Results.fSCA(k);
        end
        newSt.fSCA = tmpSCA/sum(tmpSCA);
    end
end
end