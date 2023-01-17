function [newSt] = fillSubstrate(Spectrum,Results)
%fill substrate structure with values from R0, at wavelengths in Spectrum,
%singles converted to doubles to accommodate fmincon
R0 = Results.r0;
if isempty(R0)
    R0 = 0;
    % function will fill below
end

% if R0 is an interpolant
if contains(class(R0),'cfit') ||...
        contains(class(R0),'interpolant','IgnoreCase',true)
    newSt.fR0 = {R0};
    
    % R0 is a cell vector of interpolants
elseif iscell(R0)
    for k=1:length(R0)
        assert(contains(class(R0{k}),'cfit') ||...
            contains(class(R0{k}),'interpolant','IgnoreCase',true),...
            'if ''R0'' is a cell vector, cells must contain interpolating or cfit objects')
    end
    newSt.fR0 = R0(:).'; % make sure it's a row vector
    
    % R0 is a table with columns for wavelength and reflectance, set an
    % interpolant for each column
elseif istable(R0)
%     assert(all(R0.wavelength>0) && all(R0.reflectance(:)>=0) &&...
%         all(R0.reflectance(:)<=1),'R0 table bad wavelength or reflectance')
    % make sure wavelength units for R0 match those for wavelength
    assert(~isempty(R0.Properties.VariableUnits),...
        'if R0 is a table, wavelength units must be specified')
    R0.wavelength = double(convertLengthUnits(R0.wavelength,...
        R0.Properties.VariableUnits{1},Spectrum.waveUnit));
    R0.Properties.VariableUnits{1} = Spectrum.waveUnit;
    % interpolate
    newSt.fR0 = cell(1,size(R0.reflectance,2));
    for k=1:size(R0.reflectance,2)
        F = griddedInterpolant(double(R0.wavelength),...
            double(R0.reflectance(:,k)),'pchip','nearest');
        newSt.fR0{k} = F;
    end
    
    % R0 a single number or a row vector
elseif (isscalar(R0) || isrow(R0)) && isnumeric(R0)
    x = [min(Spectrum.wavelength) max(Spectrum.wavelength)];
    if isscalar(R0)
        y = double([R0 R0]);
        F = fit(x',y','poly1');
        newSt.fR0{1} = F;
    else
        for k=1:length(R0)
            y = double([R0(k) R0(k)]);
            F = fit(x',y','poly1');
            newSt.fR0{k} = F;
        end
    end
else
    error('class(R0) = %s, not recognized',class(R0))
end
newSt.waveUnit = Spectrum.waveUnit;
end